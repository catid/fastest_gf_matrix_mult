#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include "Platform.hpp"
#include "BitMath.hpp"
#include "AbyssinianPRNG.hpp"
#include "Clock.hpp"
using namespace cat;

static Clock m_clock;

static const int GF_BITS = 5;
static const int GF_SIZE = 1 << GF_BITS;

static u8 GEN_POLY[GF_SIZE];
static int GEN_POLY_COUNT = 0;

static void FindGeneratorPolynomials()
{
    GEN_POLY_COUNT = 0;
    cout << "static const u8 GEN_POLY[] = {" << endl;
    int seen = 0;
    for (u32 taps = 0; taps < GF_SIZE; ++taps)
    {
        u32 lfsr = 1;

        int count = 0;
        for (int ii = 0; ii < GF_SIZE - 1; ++ii)
        {
            u32 lsb = lfsr & 1;
            lfsr >>= 1;
            if (lsb) lfsr ^= taps;
            if (lfsr == 1) ++count;
        }

        if (lfsr == 1 && count == 1)
        {
            GEN_POLY[GEN_POLY_COUNT++] = taps;
            cout << "0x" << hex << (int)taps << dec << ", ";
            if ((++seen & 7) == 0) cout << endl;
        }
    }
    cout << "};" << endl;
}

static u16 GF_LOG_TABLE[GF_SIZE];
static u8 GF_EXP_TABLE[GF_SIZE * 2 * 2 + 1];
static u8 GF_INV_TABLE[GF_SIZE];

u8 * CAT_RESTRICT GF_MUL_TABLE = 0;
u8 * CAT_RESTRICT GF_DIV_TABLE = 0;

static u8 CAUCHY_ONES[GF_SIZE];

static u8 MINWEIGHT_TABLE[GF_SIZE];

// return x * y in GF
// For repeated multiplication by a constant, it is faster to put the constant in y.
static CAT_INLINE u8 GFMultiply(u8 x, u8 y) {
	return GF_MUL_TABLE[((u32)y << GF_BITS) + x];
}

// return x / y in GF
// Memory-access optimized for constant divisors in y.
static CAT_INLINE u8 GFDivide(u8 x, u8 y) {
    return GF_DIV_TABLE[((u32)y << GF_BITS) + x];
}

void InitInvTable() {
    for (int x = 0; x < GF_SIZE; ++x)
    {
		GF_INV_TABLE[x] = GFDivide(1, x);
	}
}

// Unpack multiplication tables
void InitMulDivTables() {
	// If not initialized already,
    if (GF_MUL_TABLE)
    {
        free(GF_MUL_TABLE);
    }
    GF_MUL_TABLE = (u8 *)malloc(GF_SIZE * GF_SIZE * 2);
    GF_DIV_TABLE = GF_MUL_TABLE + GF_SIZE * GF_SIZE;

	// Allocate table memory < 65KB x 2
	u8 *m = GF_MUL_TABLE, *d = GF_DIV_TABLE;

	// Unroll y = 0 subtable
    for (int x = 0; x < GF_SIZE; ++x)
    {
		m[x] = d[x] = 0;
	}

	// For each other y value,
    for (int y = 1; y < GF_SIZE; ++y)
    {
		// Calculate log(y) for mult and GF_SIZE - 1 - log(y) for div
		const u8 log_y = (u8)GF_LOG_TABLE[y];
        const u8 log_yn = GF_SIZE - 1 - log_y;

		// Next subtable
        m += GF_SIZE;
        d += GF_SIZE;

		// Unroll x = 0
		m[0] = 0;
		d[0] = 0;

		// Calculate x * y, x / y
        for (int x = 1; x < GF_SIZE; ++x)
        {
			int log_x = GF_LOG_TABLE[x];

			m[x] = GF_EXP_TABLE[log_x + log_y];
			d[x] = GF_EXP_TABLE[log_x + log_yn];
		}
	}
}

void GenerateExpLogTables(int ii, u16 LogTable[GF_SIZE], u8 ALogTable[GF_SIZE * 2 * 2 + 1])
{
	u32 poly = (GEN_POLY[ii] << 1) | 1;

    LogTable[0] = GF_SIZE * 2;
	ALogTable[0] = 1;
    for (u32 jj = 1; jj < GF_SIZE - 1; ++jj)
	{
		u32 next = (u32)ALogTable[jj - 1] * 2;
        if (next >= GF_SIZE) next ^= poly;

		ALogTable[jj] = next;
		LogTable[ALogTable[jj]] = jj;
	}

    ALogTable[GF_SIZE - 1] = ALogTable[0];
    LogTable[ALogTable[GF_SIZE - 1]] = GF_SIZE - 1;

    for (u32 jj = GF_SIZE; jj < 2 * (GF_SIZE - 1); ++jj)
	{
        ALogTable[jj] = ALogTable[jj % (GF_SIZE - 1)];
	}

    ALogTable[2 * (GF_SIZE - 1)] = 1;

    for (u32 jj = 2 * (GF_SIZE - 1) + 1; jj < 4 * (GF_SIZE - 1); ++jj)
	{
		ALogTable[jj] = 0;
	}
}

// Calculate number of 1s in Cauchy 8x8 submatrix representation
static int cauchy_ones(u8 n) {
	/*
	 * w = 8 so the Cauchy representation is an 8x8 submatrix
	 * in place of the GF(256) values of the matrix.
	 *
	 * To generate the 8x8 submatrix, the first column is the
	 * original value in binary.  And the remaining columns
	 * are the column to the left times 2 in GF(256).
	 */

	// Count number of ones in the first column
	int ones = BIT_COUNT_TABLE[n];

	// For each remaining column,
	for (int w = 1; w < GF_BITS; ++w) {
		// Generate the column
		// NOTE: Optimizes to a single 256-byte table lookup
		n = GFMultiply(n, 2);

		// Count the bits in it
		ones += BIT_COUNT_TABLE[n];
	}

	return ones;
}

static void bitrow_print(u8 n)
{
    u8 mask = 1;
    for (int i = 0; i < GF_BITS; ++i)
    {
        if (n & mask)
        {
            cout << "1";
        }
        else
        {
            cout << "0";
        }

        mask <<= 1;
    }
    cout << endl;
}

static void bitmatrix_print(u8 n)
{
    cout << endl << "Bit pattern for " << hex << (int)n << dec << endl;
    bitrow_print(n);

    // For each remaining column,
    for (int w = 1; w < GF_BITS; ++w)
    {
        // Generate the column
        // NOTE: Optimizes to a single 256-byte table lookup
        n = GFMultiply(n, 2);

        bitrow_print(n);
    }
}

void GenerateCauchyOnes() {
    for (int x = 0; x < GF_SIZE; ++x)
    {
		CAUCHY_ONES[x] = cauchy_ones(x);
	}
}

void print(int k, int m, u8 *matrix, bool cstyle = false) {
	if (cstyle) {
		cout << "static const u8 CAUCHY_MATRIX_" << m << "[" << (m - 1) << " * " << k << "] = {" << endl;
		for (int y = m > 1 ? 1 : 0; y < m; ++y) {
			if (y > 1) {
				cout << "// For row " << y << ":" << endl;
			}
			for (int x = 0; x < k; ++x) {
				cout << dec << (int)matrix[y * k + x];
				if (x > 0 && (x % 20) == 19) {
					cout << "," << endl;
				} else if (x < k - 1) {
					cout << ",";
				}
			}
			if (y < m - 1) {
				cout << "," << endl;
			}
		}
		cout << dec << "};" << endl;
	} else {
		cout << "[" << endl;
		for (int y = 0; y < m; ++y) {
			for (int x = 0; x < k; ++x) {
				cout << hex << setw(2) << setfill('0') << (int)matrix[y * k + x] << " ";
			}
			cout << endl;
		}
		cout << dec << "]" << endl;
    
        cout << "{" << endl;
        for (int y = 0; y < m; ++y)
        {
            for (int w = 0; w < GF_BITS; ++w)
            {
                for (int x = 0; x < k; ++x)
                {
                    u8 n = matrix[y * k + x];

                    for (int i = 0; i < w; ++i)
                    {
                        n = GFMultiply(n, 2);
                    }

                    u8 mask = 1;
                    for (int i = 0; i < GF_BITS; ++i)
                    {
                        if (n & mask)
                        {
                            cout << "1";
                        }
                        else
                        {
                            cout << "0";
                        }
                        mask <<= 1;
                    }
                    cout << " ";
                }

                cout << endl;
            }

            cout << endl;
        }
        cout << dec << "}" << endl;
    }
}

void SortColumns(int k, int m, u8 *matrix) {
	int *counts = new int[k];

	for (int x = 0; x < k; ++x) {
		int ones = 0;
		for (int y = 0; y < m; ++y) {
			ones += CAUCHY_ONES[matrix[y*k + x]];
		}
		counts[x] = ones;
	}

	for (int x = 0; x < k; ++x) {
		int smallest = counts[x], best_x = x;

		for (int z = x + 1; z < k; ++z) {
			int ones = counts[z];

			if (counts[z] < smallest) {
				smallest = ones;
				best_x = z;
			}
		}

		// swap counts
		counts[best_x] = counts[x];
		counts[x] = smallest;

		// swap columns
		for (int y = 1; y < m; ++y) {
			u8 temp = matrix[y*k + x];
			matrix[y*k + x] = matrix[y*k + best_x];
			matrix[y*k + best_x] = temp;
		}
	}

	delete []counts;
}

/*
 * Cauchy matrices are defined by two vectors X, Y s.t. X, Y share no elements
 * in common from the set GF(256).
 *
 * Each element i,j of the Cauchy matrix is 1/(Xi + Yj).
 *
 * Another useful property is that you can multiply each row or column of a
 * Cauchy matrix and it will still be invertible.
 *
 * Since the number of Cauchy ones is far better for 1 than other elements of
 * GF(256), all of the best options will have a 1 in at least one row of each
 * column.  And one column can be all 1s.
 */

int GenerateCauchyMatrix(int k, int m, u8 *matrix, u8 *X, u8 *Y) {
	int ones = CAUCHY_ONES[1] * k;

	u8 *row = matrix;

	for (int y = 0; y < k; ++y) {
		row[y] = 1;
	}

	for (int y = 1; y < m; ++y) {
		u8 YC = Y[y];

		row += k;

		for (int x = 0; x < k; ++x) {
			u8 XC = X[x];

			u8 D = Y[0] ^ XC;
			u8 C = GFMultiply(GF_INV_TABLE[XC ^ YC], D);

			row[x] = C;
		}
	}

	row = matrix;
	for (int y = 1; y < m; ++y) {
		row += k;

		int best = 0x7fffffff, best_x = 0;

		for (int x = 0; x < k; ++x) {
			u8 XC = row[x];

			int count = 0;

			for (int z = 0; z < k; ++z) {
				u8 C = GFDivide(row[z], XC);

				count += CAUCHY_ONES[C];
			}

			if (count < best) {
				best = count;
				best_x = x;
			}
		}

		u8 XC = row[best_x];

		for (int z = 0; z < k; ++z) {
			u8 C = GFDivide(row[z], XC);

			row[z] = C;

			ones += CAUCHY_ONES[C];
		}
	}

	return ones;
}

int ImproveMatrixRows(int k, int subk, int m, u8 *matrix) {
	for (int y = 1; y < m; ++y) {
		int best = 0x7fffffff, best_A;
		for (int x = 0; x < k; ++x) {
			u8 A = matrix[y*k + x];
			u8 IA = GF_INV_TABLE[A];

			int ones = 0;
			for (int z = 0; z < k; ++z) {
				u8 B = matrix[y*k + z];
				u8 M = GFMultiply(B, IA);
				ones += CAUCHY_ONES[M];
			}
			if (ones < best) {
				best = ones;
				best_A = IA;
			}
		}
		for (int z = 0; z < k; ++z) {
			u8 B = matrix[y*k + z];
			u8 M = GFMultiply(B, best_A);
			matrix[y*k + z] = M;
		}
	}

	SortColumns(k, m, matrix);

	int total = 0;
	for (int y = 0; y < m; ++y) {
		for (int x = 0; x < subk; ++x) {
			total += CAUCHY_ONES[matrix[y*k + x]];
		}
	}
	return total;
}

void SolveBestMatrix(int m, int subk) {
	//   A B C D E
	// F 1 1 1 1 1
	// G a b c d e
	// H f g h i j

    const int k = GF_SIZE - m;
	u8 *matrix = new u8[k * m];
	u8 *best_matrix = new u8[k * m];
	int best_matrix_ones = 0x7fffffff;
    u8 X[GF_SIZE], Y[GF_SIZE];

	double t0 = m_clock.usec();

	// First row is always all ones
	for (int x = 0; x < k; ++x) {
		matrix[x] = 1;
	}

	// Choose a seed of A,F and solve the rest with a greedy algorithm
	//for (int F = 0; F < 256; ++F) {
		//for (int A = 0; A < 256; ++A) {
			//if (A != F) {
	int F = 0, A = 1; {{{
				u8 seen[GF_SIZE];
                for (int ii = 0; ii < GF_SIZE; ++ii)
                {
					seen[ii] = 0;
				}
				seen[A] = 1;
				seen[F] = 1;
				X[0] = A;
				Y[0] = F;
				u8 AF = A ^ F;

				// In the order of greedy solution solve the Y values first.

                int trial_ones = GF_BITS * subk;

				// Unroll the first column
				for (int y = 1; y < m; ++y) {
					// Pick the next worst element of GF(256) in weight
                    for (int ii = 1; ii < GF_SIZE; ++ii)
                    {
						u8 a = MINWEIGHT_TABLE[ii];

						// a * (A + G) = A + F
						// aA + aG = A + F
						// G = (A + F + aA) / a
						u8 G = GFDivide(AF ^ GFMultiply(a, A), a);

						if (seen[G]) {
							continue;
						}

						seen[G] = 1;
						Y[y] = G;
						matrix[y * k] = a;
						trial_ones += CAUCHY_ONES[a];
						break;
					}
				}

				// Now solve the X values.
				// TODO: Maybe iterate through this inside Y value solutions

				// For each remaining column,
				for (int x = 1; x < k; ++x) {
					int best_ones = 0x7fffffff, best_B;

					// Verify that a solution is possible for all column values
                    for (int B = 0; B < GF_SIZE; ++B)
                    {
						if (seen[B]) {
							continue;
						}

						int ones = 0;
						for (int y = 1; y < m; ++y) {
							u8 b = GFDivide(B ^ F, B ^ Y[y]);
							ones += CAUCHY_ONES[b];
						}
						if (ones < best_ones) {
							best_ones = ones;
							best_B = B;
						}
					}

					int B = best_B;
					X[x] = B;
					seen[B] = 1;

					for (int y = 1; y < m; ++y) {
						u8 b = GFDivide(B ^ F, B ^ Y[y]);
						matrix[y * k + x] = b;
					}
					if (x < subk) {
						trial_ones += best_ones;
					}
				}
				//print(k, m, matrix);

				double t1 = m_clock.usec();

				int improved_ones = ImproveMatrixRows(k, subk, m, matrix);

				double t2 = m_clock.usec();

				if (improved_ones < best_matrix_ones) {
					best_matrix_ones = improved_ones;
					memcpy(best_matrix, matrix, k * m);

					cout << "Pre-improved ones = " << trial_ones << " in " << (t1 - t0) << " usec" << endl;
					cout << "Best ones for first " << subk << " columns = " << best_matrix_ones << " in " << (t2 - t0) << " usec" << endl;
					print(k, m, best_matrix, true);

					return;
				}
			}
		}
	}
}

void SortMinWeightElements(u8 *elements) {
    int *counts = new int[GF_SIZE];

    for (int x = 0; x < GF_SIZE; ++x)
    {
		counts[x] = CAUCHY_ONES[x];
	}

    for (int x = 0; x < GF_SIZE; ++x)
    {
		int smallest = counts[x], best_x = x;

        for (int z = x + 1; z < GF_SIZE; ++z)
        {
			int ones = counts[z];

			if (counts[z] < smallest) {
				smallest = ones;
				best_x = z;
			}
		}

		// swap counts
		counts[best_x] = counts[x];
		counts[x] = smallest;

		// swap columns
		u8 temp = elements[best_x];
		elements[best_x] = elements[x];
		elements[x] = temp;
	}

	delete []counts;
}

static void ShuffleDeck8(Abyssinian &prng, u8 * CAT_RESTRICT deck)
{
	deck[0] = 0;
	const int count = GF_SIZE;

	for (u32 ii = 1;;)
	{
		u32 jj, rv = prng.Next();

		// 8-bit unroll
		switch (count - ii)
		{
			default:
				jj = (u8)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 8) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 16) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 24) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				break;

			case 3:
				jj = (u8)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
			case 2:
				jj = (u8)(rv >> 8) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
			case 1:
				jj = (u8)(rv >> 16) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
			case 0:
				return;
		}
	}
}

int PrintMinWeights() {
	const int k = GF_SIZE;
	const int m = 2;
	u8 *matrix = new u8[k * m];

    int bestPoly = 0;
    int bestCount = 32767;

	for (int ii = 0; ii < GEN_POLY_COUNT; ++ii) {
		cout << "*** For generator " << ii << ":" << endl;

		GenerateExpLogTables(ii, GF_LOG_TABLE, GF_EXP_TABLE);

		InitMulDivTables();

		InitInvTable();

		GenerateCauchyOnes();

		for (int x = 0; x < GF_SIZE; ++x) {
			matrix[x] = 1;
            matrix[x + GF_SIZE] = x;
		}

		cout << "Symbols in order:" << endl;

		SortColumns(k, m, matrix);

		print(k, m, matrix, true);

        for (int x = 1; x <= 32 && x <= GF_SIZE; ++x)
        {
			int ones = 0;
			for (int z = 1; z <= x; ++z) {
				ones += CAUCHY_ONES[matrix[k + z]];
			}

            if (x == 4 && ones < bestCount)
            {
                bestCount = ones;
                bestPoly = ii;
                cout << "-------------BEST POLY SO FAR" << endl;
            }

			cout << x << " columns = " << ones << " ones" << endl;
		}

//         for (int n = 0; n < GF_SIZE; ++n)
//         {
//             bitmatrix_print((u8)n);
//         }
	}

    return bestPoly;
}

int count_xors(int k, int m, u8* matrix, int column)
{
    bool seen_pattern[GF_SIZE] = {};
    seen_pattern[0] = true; // 0
    for (int i = 0; i < GF_BITS; ++i)
    {
        seen_pattern[1 << i] = true; // powers of 2 are free
    }

    int pattern_count = 0;

    for (int r = 0; r < m; ++r)
    {
        const int c = column;
        {
            // grab value
            u8 n = matrix[r * k + c];

            if (!seen_pattern[n])
            {
                seen_pattern[n] = true;
                ++pattern_count;
            }

            // For each remaining column,
            for (int w = 1; w < GF_BITS; ++w)
            {
                // Generate the column
                // NOTE: Optimizes to a single 256-byte table lookup
                n = GFMultiply(n, 2);

                if (!seen_pattern[n])
                {
                    seen_pattern[n] = true;
                    ++pattern_count;
                }
            }
        }
    }

    // Now find minimum number of intermediate additional XORs to do to branch between symbols

    static vector<u8> counts[GF_BITS + 1];
    int maxCount = 0;

    for (int i = 0; i <= GF_BITS; ++i)
    {
        counts[i].clear();
    }

    for (int i = 0; i < GF_SIZE; ++i)
    {
        if (seen_pattern[i])
        {
            int count = BitCount<u8>(i);

            counts[count].push_back((u8)i);

            if (count > maxCount)
            {
                maxCount = count;
            }
        }
    }

    // For counts at or above 3 we need to chart a course
    for (int i = 3; i <= maxCount; ++i)
    {
        int size = (int)counts[i].size();
        for (int j = 0; j < size; ++j)
        {
            // what we want to build
            u8 target = counts[i][j];

            // find the closest other symbol in hamming distance
            int best_diff = 1000;
            u8 best_cand = 0;

            // from all other symbols with fewer bits set
            for (int q = i - 1; q >= 2; --q)
            {
                int psize = (int)counts[q].size();
                for (int p = 0; p < psize; ++p)
                {
                    u8 cand = counts[q][p];

                    int diff = BitCount<u8>(cand ^ target);

                    if (diff <= best_diff)
                    {
                        best_diff = diff;
                        best_cand = cand;
                    }
                }
            }

            // If we found a direct path,
            if (best_cand)
            {
                if (best_diff <= 1)
                {
                    // great!  this means that we only really needed one xor to get here,
                    // and the pattern count already reflects the number of xors needed.
                }
                else
                {
                    // GREEDY HEURISTIC WARNING
                    // not so great.  we needed to xor more than once.

                    // I'm currently not xoring symbols together like 1100 + 0011 to get 1111
                    // so it would count 2 instead of 1 in this case and overestimate complexity
                    // incorrectly.  I'm hoping that this won't make a huge difference in real
                    // performance of the generated chain.  It would be worthwhile to improve
                    // the heuristic and try all these options.

                    // My goal is to have the multiplication code follow a set of predetermined
                    // instructions for constructing each row, so maybe if a heuristic is used
                    // there it should go here too.

                    // In this case, we should add intermediate patterns that do not exist
                    // towards the target.  This will help future searches also.
                    pattern_count += best_diff - 1;
                }
            }
            else
            {
                // APPROXIMATION WARNING
                // we don't have anything that remotely helps reduce the pattern count.
                // In this case, we should add intermediate patterns that do not exist
                // towards the target.  This will help future searches also.
                pattern_count += BitCount<u8>(target) - 1;
            }
        }
    }

    return pattern_count;
}

void Explore(int k, int m) {
	u8 *matrix = new u8[k * m];
	u8 *best = new u8[k * m];
    u8 *XY = new u8[GF_SIZE];
	u8 *X = XY;
	u8 *Y = XY + k;

	Abyssinian prng;

	prng.Initialize(m_clock.cycles());

	int least = 32767;

	for (int ii = 0; ii < 1000000000; ++ii) {
		ShuffleDeck8(prng, XY);

		GenerateCauchyMatrix(k, m, matrix, X, Y);

        // number of temporary variables to keep
        int pattern_count = 0;
        for (int j = 0; j < k; ++j)
        {
            pattern_count += count_xors(k, m, matrix, j);
        }

        if (pattern_count < least)
        {
            least = pattern_count;
			memcpy(best, matrix, k * m);

            // normalized cost = number of additional xors
            float cost = pattern_count / (float)GF_BITS;

            // this is the count of fractional xors in excess of a pure parity scheme.
            // the ideal would be just xoring all the inputs together, but that is only
            // achieved by the first matrix row.  an ideal matrix multiplier would reuse
            // any temporary xors generated during computation to avoid doing those xors
            // more than once.  if the temporary xor count is low then this is feasible.

            cout << "Found a better matrix with T = " << pattern_count << " temp variables and C = " << cost << " normalized cost:" << endl;
			print(k, m, best);
		}
	}

	SortColumns(k, m, best);

	cout << "Sorted matrix:" << endl;

	print(k, m, best);

	delete []matrix;
	delete []best;
	delete []XY;
}

void PrintTables() {
    print(GF_SIZE *2* 2 + 1, 1, GF_EXP_TABLE, true);
    print(GF_SIZE, 1, GF_INV_TABLE, true);

	cout << "static const u8 GFC_LOG[] = {" << endl;
	for (int x = 0; x < GF_SIZE; ++x) {
		cout << dec << (int)GF_LOG_TABLE[x];
		if (x > 0 && (x % 20) == 19) {
			cout << "," << endl;
		} else {
			cout << ",";
		}
	}
	cout << dec << "};" << endl;
}

int main() {
	m_clock.OnInitialize();

	cout << "Exploring options..." << endl;

    FindGeneratorPolynomials();

    int bestPoly = PrintMinWeights();

    GenerateExpLogTables(bestPoly, GF_LOG_TABLE, GF_EXP_TABLE);

	InitMulDivTables();

	InitInvTable();

	GenerateCauchyOnes();

	for (int ii = 0; ii < GF_SIZE; ++ii) {
		MINWEIGHT_TABLE[ii] = ii;
	}
	SortMinWeightElements(MINWEIGHT_TABLE);

    cout << "Min Weight table:" << endl;
    print(GF_SIZE, 1, MINWEIGHT_TABLE, true);

	//SolveBestMatrix(4, 4);
	Explore(4, 4);

	PrintTables();

	m_clock.OnFinalize();

	return 0;
}

