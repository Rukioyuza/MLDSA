--!native
--!optimize 2
-- Polynomials in Z_q[X]/(X²⁵⁶ + 1), q = 8380417 (FIPS 204 §7.5). Immutable: every operation returns a new polynomial.
-- A polynomial is a table of 256 coefficients in [0, q), indexed 0..255 as in the specification.
-- `+`, `-` and unary `-` work in either domain; `*` multiplies two polynomials that are both in the NTT domain.
-- Products of two coefficients stay below 2⁴⁶, so plain doubles are exact.

local MODULUS = 8380417 -- q
local ROOT = 1753 -- ζ, a primitive 512th root of unity mod q

local Polynomial = {}
Polynomial.__index = Polynomial

local function build(coefficient: (index: number) -> number) -- Σ coefficient(i)·Xⁱ, each reduced mod q
	local polynomial = setmetatable({}, Polynomial)
	for index = 0, 255 do polynomial[index] = coefficient(index) % MODULUS end
	return polynomial
end
Polynomial.new = build

local function centered(value: number, modulus: number): number -- value mod± modulus, in (−modulus/2, modulus/2]
	local reduced = value % modulus
	return if reduced > modulus // 2 then reduced - modulus else reduced
end

---- Arithmetic ----

function Polynomial.__add(left, right)
	return build(function(index) return left[index] + right[index] end)
end

function Polynomial.__sub(left, right)
	return build(function(index) return left[index] - right[index] end)
end

function Polynomial.__unm(polynomial)
	return build(function(index) return -polynomial[index] end)
end

function Polynomial.__mul(left, right) -- MultiplyNTT (Algorithm 45): pointwise
	return build(function(index) return left[index] * right[index] end)
end

function Polynomial:norm(): number -- ‖·‖∞ of the centred coefficients
	local largest = 0
	for index = 0, 255 do largest = math.max(largest, math.abs(centered(self[index], MODULUS))) end
	return largest
end

local function bitReverse8(value: number): number
	local reversed = 0
	for _ = 1, 8 do
		reversed, value = 2 * reversed + value % 2, value // 2
	end
	return reversed
end

local zetas = {} -- zetas[i] = ζ^BitRev8(i)
do
	local powers = { [0] = 1 }
	for exponent = 1, 255 do powers[exponent] = powers[exponent - 1] * ROOT % MODULUS end
	for index = 0, 255 do zetas[index] = powers[bitReverse8(index)] end
end

function Polynomial:ntt() -- NTT (Algorithm 41)
	local result, zetaIndex, length = table.clone(self), 1, 128
	while length >= 1 do
		for start = 0, 255, 2 * length do
			local zeta = zetas[zetaIndex]
			zetaIndex += 1
			for index = start, start + length - 1 do
				local low, high = result[index], zeta * result[index + length]
				result[index], result[index + length] = (low + high) % MODULUS, (low - high) % MODULUS
			end
		end
		length //= 2
	end
	return result
end

function Polynomial:inverseNtt() -- NTT⁻¹ (Algorithm 42); 8347681 = 256⁻¹ mod q
	local result, zetaIndex, length = table.clone(self), 255, 1
	while length < 256 do
		for start = 0, 255, 2 * length do
			local zeta = zetas[zetaIndex]
			zetaIndex -= 1
			for index = start, start + length - 1 do
				local low, high = result[index], result[index + length]
				result[index], result[index + length] = (low + high) % MODULUS, zeta * (high - low) % MODULUS
			end
		end
		length *= 2
	end
	return build(function(index) return result[index] * 8347681 end)
end

---- Rounding (§7.4) ----

function Polynomial:power2Round() -- Power2Round (Algorithm 35): t = t₁·2¹³ + t₀, t₀ ∈ (−2¹², 2¹²]
	return build(function(index) return (self[index] - centered(self[index], 8192)) // 8192 end),
		build(function(index) return centered(self[index], 8192) end)
end

local function decompose(value: number, gamma2: number): (number, number) -- Decompose (Algorithm 36): r = r₁·2γ₂ + r₀
	value %= MODULUS
	local low = centered(value, 2 * gamma2)
	if value - low == MODULUS - 1 then return 0, low - 1 end
	return (value - low) // (2 * gamma2), low
end

function Polynomial:highBits(gamma2: number) -- HighBits (Algorithm 37)
	return build(function(index) return (decompose(self[index], gamma2)) end)
end

function Polynomial:lowBits(gamma2: number) -- LowBits (Algorithm 38)
	return build(function(index) return select(2, decompose(self[index], gamma2)) end)
end

function Polynomial.makeHint(z, r, gamma2: number) -- MakeHint (Algorithm 39): 1 where adding z changes the high bits of r
	return build(function(index)
		return if decompose(r[index], gamma2) == decompose(r[index] + z[index], gamma2) then 0 else 1 end)
end

function Polynomial.useHint(hint, r, gamma2: number) -- UseHint (Algorithm 40): high bits of r, corrected by the hint
	local count = (MODULUS - 1) // (2 * gamma2)
	return build(function(index)
		local high, low = decompose(r[index], gamma2)
		if hint[index] == 0 then return high end
		return (if low > 0 then high + 1 else high - 1) % count
	end)
end

---- Encoding (§7.1): each coefficient c as `bits` bits, packed little-endian ----
---- SimpleBitPack stores c itself; BitPack (when `upper` is given) stores upper − c ----

function Polynomial:encode(bits: number, upper: number?): buffer -- → 32·bits bytes
	local bytes = buffer.create(32 * bits)
	for index = 0, 255 do
		buffer.writebits(bytes, index * bits, bits, if upper then (upper - self[index]) % MODULUS else self[index])
	end
	return bytes
end

function Polynomial.decode(bytes: buffer, offset: number, bits: number, upper: number?) -- SimpleBitUnpack / BitUnpack
	return build(function(index)
		local value = buffer.readbits(bytes, 8 * offset + index * bits, bits)
		return if upper then upper - value else value
	end)
end

---- Sampling ----

function Polynomial.sampleUniform(stream) -- RejNTTPoly (Algorithm 30): keep 23-bit candidates below q
	local accepted = {}
	while #accepted < 256 do
		local candidate = buffer.readbits(stream:squeeze(3), 0, 23)
		if candidate < MODULUS then table.insert(accepted, candidate) end
	end
	return build(function(index) return accepted[index + 1] end)
end

function Polynomial.sampleBounded(stream, eta: number) -- RejBoundedPoly (Algorithm 31): half-bytes mapped into [−η, η]
	local accepted = {}
	while #accepted < 256 do
		local byte = buffer.readu8(stream:squeeze(1), 0)
		for _, half in { byte % 16, byte // 16 } do -- CoeffFromHalfByte (Algorithm 15)
			if #accepted == 256 then break end
			if eta == 2 and half < 15 then
				table.insert(accepted, 2 - half % 5)
			elseif eta == 4 and half < 9 then
				table.insert(accepted, 4 - half)
			end
		end
	end
	return build(function(index) return accepted[index + 1] end)
end

function Polynomial.sampleInBall(stream, tau: number) -- SampleInBall (Algorithm 29): τ coefficients ±1, the rest 0
	local signs, coefficients = stream:squeeze(8), {}
	for index = 0, 255 do coefficients[index] = 0 end
	for index = 256 - tau, 255 do
		local swap
		repeat swap = buffer.readu8(stream:squeeze(1), 0) until swap <= index
		coefficients[index] = coefficients[swap]
		coefficients[swap] = 1 - 2 * buffer.readbits(signs, index + tau - 256, 1)
	end
	return build(function(index) return coefficients[index] end)
end

return Polynomial
