--!native
--!optimize 2
-- ML-DSA (FIPS 204), deterministic core: the caller supplies all randomness as 32-byte buffers.
-- Signing is pure ML-DSA with an optional context string (≤ 255 bytes); all-zero randomness gives the deterministic variant.
local Keccak = require("@self/Keccak")
local Polynomial = require("@self/Polynomial")

---- Bytes ----

local function slice(source: buffer, offset: number, length: number): buffer
	local part = buffer.create(length)
	buffer.copy(part, 0, source, offset, length)
	return part
end

local function concat(...: buffer): buffer
	local length, offset = 0, 0
	for _, part in { ... } do length += buffer.len(part) end
	local joined = buffer.create(length)
	for _, part in { ... } do
		buffer.copy(joined, offset, part)
		offset += buffer.len(part)
	end
	return joined
end

local function equal(left: buffer, right: buffer): boolean -- same lengths; inspects every byte, no early exit
	local difference = 0
	for index = 0, buffer.len(left) - 1 do
		difference = bit32.bor(difference, bit32.bxor(buffer.readu8(left, index), buffer.readu8(right, index)))
	end
	return difference == 0
end

local function bitLength(value: number): number
	local length = 0
	while value > 0 do length, value = length + 1, value // 2 end
	return length
end

---- Hash functions (FIPS 204 §3.7): arguments are absorbed in order, each a buffer or a single byte ----

local H, G = Keccak.shake256, Keccak.shake128 -- both return the sponge; call :squeeze(n)

local function messageHash(tr, message, context) -- μ = H(tr ‖ M′, 64) with M′ = 0 ‖ |ctx| ‖ ctx ‖ M (Algorithms 2-3)
	return H(tr, 0, buffer.len(context), context, message):squeeze(64)
end

---- Vectors and matrices of polynomials (plain arrays) ----

local function vectorOf(rank: number, entry: (index: number) -> any)
	local entries = table.create(rank)
	for index = 1, rank do entries[index] = entry(index) end
	return entries
end

local function dot(left, right) -- Σ left[i] ∘ right[i], both in the NTT domain
	local sum = left[1] * right[1]
	for index = 2, #left do sum += left[index] * right[index] end
	return sum
end

local function ntt(vector)
	return vectorOf(#vector, function(index) return vector[index]:ntt() end)
end

local function norm(vector): number -- ‖·‖∞ over every coefficient
	local largest = 0
	for _, polynomial in vector do largest = math.max(largest, polynomial:norm()) end
	return largest
end

local function encodeVector(vector, bits: number, upper: number?): buffer
	return concat(table.unpack(vectorOf(#vector, function(index) return vector[index]:encode(bits, upper) end)))
end

local function decodeVector(rank: number, bytes: buffer, offset: number, bits: number, upper: number?)
	return vectorOf(rank, function(index)
		return Polynomial.decode(bytes, offset + (index - 1) * 32 * bits, bits, upper)
	end)
end

---- Scheme ----

local function parameterSet(rows, columns, eta, tau, lambda, gamma1, gamma2, omega) -- (k, ℓ, η, τ, λ, γ₁, γ₂, ω)
	local beta, challengeBytes = tau * eta, lambda // 4
	local etaBits, maskBits = bitLength(2 * eta), 1 + bitLength(gamma1 - 1)
	local highBits = bitLength((8380417 - 1) // (2 * gamma2) - 1) -- bits per coefficient of w₁
	local publicKeyBytes = 32 + 320 * rows
	local s2Offset = 128 + 32 * etaBits * columns
	local t0Offset = s2Offset + 32 * etaBits * rows
	local secretKeyBytes = t0Offset + 32 * 13 * rows
	local hintOffset = challengeBytes + 32 * maskBits * columns
	local signatureBytes = hintOffset + omega + rows

	local function expandMatrix(rho) -- ExpandA (Algorithm 32): Â[r][s] = RejNTTPoly(ρ ‖ s ‖ r), with r and s counted from 0
		return vectorOf(rows, function(row)
			return vectorOf(columns, function(column)
				return Polynomial.sampleUniform(G(rho, column - 1, row - 1))
			end)
		end)
	end

	local function expandSecrets(seed) -- ExpandS (Algorithm 33): s₁ from counters 0..ℓ−1, s₂ from ℓ..ℓ+k−1 (2 bytes each)
		local function sample(counter) return Polynomial.sampleBounded(H(seed, counter, 0), eta) end
		return vectorOf(columns, function(index) return sample(index - 1) end),
			vectorOf(rows, function(index) return sample(columns + index - 1) end)
	end

	local function expandMask(seed, kappa) -- ExpandMask (Algorithm 34): y[r] = BitUnpack(H(ρ″ ‖ κ + r, 32c), γ₁ − 1, γ₁)
		return vectorOf(columns, function(index)
			local counter = kappa + index - 1
			local bytes = H(seed, counter % 256, counter // 256 % 256):squeeze(32 * maskBits)
			return Polynomial.decode(bytes, 0, maskBits, gamma1)
		end)
	end

	local function challengeHash(mu, w1) -- c̃ = H(μ ‖ w1Encode(w₁), λ/4)
		return H(mu, encodeVector(w1, highBits)):squeeze(challengeBytes)
	end

	local function encodeHint(hint) -- HintBitPack (Algorithm 20): positions of the ones, then each polynomial's running count
		local bytes, count = buffer.create(omega + rows), 0
		for row = 1, rows do
			for index = 0, 255 do
				if hint[row][index] == 1 then
					buffer.writeu8(bytes, count, index)
					count += 1
				end
			end
			buffer.writeu8(bytes, omega + row - 1, count)
		end
		return bytes
	end

	local function decodeHint(bytes, offset) -- HintBitUnpack (Algorithm 21); nil unless canonical
		local hint, count = table.create(rows), 0
		for row = 1, rows do
			local limit = buffer.readu8(bytes, offset + omega + row - 1)
			if limit < count or limit > omega then return nil end
			local ones = {}
			for position = count, limit - 1 do
				local index = buffer.readu8(bytes, offset + position)
				if position > count and buffer.readu8(bytes, offset + position - 1) >= index then return nil end
				ones[index] = true
			end
			hint[row], count = Polynomial.new(function(index) return if ones[index] then 1 else 0 end), limit
		end
		for position = count, omega - 1 do
			if buffer.readu8(bytes, offset + position) ~= 0 then return nil end
		end
		return hint
	end

	-- ML-DSA (FIPS 204 §5-6)

	local function generateKeys(seed) -- ML-DSA.KeyGen_internal (Algorithm 6) → publicKey = ρ ‖ t₁, secretKey = ρ ‖ K ‖ tr ‖ s₁ ‖ s₂ ‖ t₀
		assert(buffer.len(seed) == 32, "bad seed length")
		local expanded = H(seed, rows, columns):squeeze(128)
		local rho, key = slice(expanded, 0, 32), slice(expanded, 96, 32)
		local matrix = expandMatrix(rho)
		local s1, s2 = expandSecrets(slice(expanded, 32, 64))
		local s1Hat, t1, t0 = ntt(s1), {}, {}
		for row = 1, rows do -- t = NTT⁻¹(Â∘ŝ₁) + s₂ = t₁·2¹³ + t₀
			t1[row], t0[row] = (dot(matrix[row], s1Hat):inverseNtt() + s2[row]):power2Round()
		end
		local publicKey = concat(rho, encodeVector(t1, 10))
		local tr = H(publicKey):squeeze(64)
		return publicKey,
			concat(rho, key, tr, encodeVector(s1, etaBits, eta), encodeVector(s2, etaBits, eta), encodeVector(t0, 13, 4096))
	end

	local function sign(secretKey, message, randomness, context: buffer?) -- ML-DSA.Sign (Algorithm 2) → signature = c̃ ‖ z ‖ h
		context = context or buffer.create(0)
		assert(buffer.len(secretKey) == secretKeyBytes and buffer.len(randomness) == 32, "bad input length")
		assert(buffer.len(context) <= 255, "context too long")
		-- ML-DSA.Sign_internal (Algorithm 7)
		local key, tr = slice(secretKey, 32, 32), slice(secretKey, 64, 64)
		local s1Hat = ntt(decodeVector(columns, secretKey, 128, etaBits, eta))
		local s2Hat = ntt(decodeVector(rows, secretKey, s2Offset, etaBits, eta))
		local t0Hat = ntt(decodeVector(rows, secretKey, t0Offset, 13, 4096))
		local matrix = expandMatrix(slice(secretKey, 0, 32))
		local mu = messageHash(tr, message, context)
		local maskSeed = H(key, randomness, mu):squeeze(64) -- ρ″

		local kappa = 0
		while true do
			local y = expandMask(maskSeed, kappa)
			kappa += columns
			local yHat = ntt(y)
			local w = vectorOf(rows, function(row) return dot(matrix[row], yHat):inverseNtt() end) -- w = NTT⁻¹(Â∘ŷ)
			local challenge = challengeHash(mu, vectorOf(rows, function(row) return w[row]:highBits(gamma2) end))
			local cHat = Polynomial.sampleInBall(H(challenge), tau):ntt()
			local z = vectorOf(columns, function(index) return y[index] + (cHat * s1Hat[index]):inverseNtt() end)
			local approximate = vectorOf(rows, function(row) return w[row] - (cHat * s2Hat[row]):inverseNtt() end) -- w − cs₂
			local lowBits = vectorOf(rows, function(row) return approximate[row]:lowBits(gamma2) end)
			if norm(z) < gamma1 - beta and norm(lowBits) < gamma2 - beta then
				local ct0 = vectorOf(rows, function(row) return (cHat * t0Hat[row]):inverseNtt() end)
				local hint = vectorOf(rows, function(row)
					return Polynomial.makeHint(-ct0[row], approximate[row] + ct0[row], gamma2)
				end)
				local ones = 0
				for _, polynomial in hint do
					for index = 0, 255 do ones += polynomial[index] end
				end
				if norm(ct0) < gamma2 and ones <= omega then
					return concat(challenge, encodeVector(z, maskBits, gamma1), encodeHint(hint))
				end
			end
		end
	end

	local function verify(publicKey, message, signature, context: buffer?): boolean -- ML-DSA.Verify (Algorithm 3)
		context = context or buffer.create(0)
		assert(buffer.len(publicKey) == publicKeyBytes, "bad input length")
		if buffer.len(signature) ~= signatureBytes or buffer.len(context) > 255 then return false end
		-- ML-DSA.Verify_internal (Algorithm 8)
		local challenge = slice(signature, 0, challengeBytes)
		local z = decodeVector(columns, signature, challengeBytes, maskBits, gamma1)
		local hint = decodeHint(signature, hintOffset)
		if not hint or norm(z) >= gamma1 - beta then return false end
		local matrix = expandMatrix(slice(publicKey, 0, 32))
		local t1 = decodeVector(rows, publicKey, 32, 10)
		local mu = messageHash(H(publicKey):squeeze(64), message, context)
		local cHat = Polynomial.sampleInBall(H(challenge), tau):ntt()
		local zHat = ntt(z)
		local w1 = vectorOf(rows, function(row) -- UseHint(h, NTT⁻¹(Â∘ẑ − ĉ∘NTT(t₁·2¹³)))
			local shifted = Polynomial.new(function(index) return 8192 * t1[row][index] end):ntt()
			return Polynomial.useHint(hint[row], (dot(matrix[row], zHat) - cHat * shifted):inverseNtt(), gamma2)
		end)
		return equal(challenge, challengeHash(mu, w1))
	end

	return {
		generateKeys = generateKeys,
		sign = sign,
		verify = verify,
		sizes = { publicKey = publicKeyBytes, secretKey = secretKeyBytes, signature = signatureBytes },
	}
end

return {
	[44] = parameterSet(4, 4, 2, 39, 128, 2 ^ 17, 95232, 80), -- γ₂ = (q − 1)/88
	[65] = parameterSet(6, 5, 4, 49, 192, 2 ^ 19, 261888, 55), -- γ₂ = (q − 1)/32
	[87] = parameterSet(8, 7, 2, 60, 256, 2 ^ 19, 261888, 75),
}
