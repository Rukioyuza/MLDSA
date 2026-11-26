<h1 align="left">ML-DSA · Luau</h1>
 
<p align="left">
<img src="https://img.shields.io/badge/Post--Quantum-Cryptography-181717?style=for-the-badge&logo=letsencrypt&logoColor=white" />
<img src="https://img.shields.io/badge/NIST-FIPS_203-181717?style=for-the-badge" />
<img src="https://img.shields.io/badge/Pure-Luau-181717?style=for-the-badge&logo=roblox&logoColor=white" />
</p>
A boring, somewhat readable implementation of ML-DSA (NIST FIPS 204) in pure Luau.
 
## About
 
ML-DSA is the post-quantum signature NIST standardised as FIPS 204.
 
This is a boring re-implementation:
 
It is not the first ML-DSA for Luau. [rbx-cryptography](https://github.com/daily3014/rbx-cryptography) has had ML-DSA-44, 65 and 87 since August 2025 and is considerably faster (see below). If you need ML-DSA in production, use that over this lying repo.

## Performance

| Operation | Time |
|---|---|
| ML-DSA-44 generateKeys | 3.2 ms |
| ML-DSA-44 sign | 6.2 ms |
| ML-DSA-44 verify | 3.6 ms |

### Compared with rbx-cryptography

| Operation | rbx-cryptography | Ratio |
|---|---|---|
| ML-DSA-44 generateKeys | 42 µs | 76x slower |
| ML-DSA-44 sign | 47 µs | 132x slower |
| ML-DSA-44 verify | 20 µs | 180x slower |
 
## Usage
Don't use it, use https://github.com/daily3014/rbx-cryptography its faster, actually passes all tests and more secure
 
