#include "Platform.hpp"
using namespace cat;



//// GF(256) math

// Tables generated with optimal polynomial 0x187 = 110000111b

static const u16 GFC256_LOG_TABLE[256] = {
    512, 255, 1, 99, 2, 198, 100, 106, 3, 205, 199, 188, 101, 126, 107, 42, 4, 141, 206, 78,
    200, 212, 189, 225, 102, 221, 127, 49, 108, 32, 43, 243, 5, 87, 142, 232, 207, 172, 79, 131,
    201, 217, 213, 65, 190, 148, 226, 180, 103, 39, 222, 240, 128, 177, 50, 53, 109, 69, 33, 18,
    44, 13, 244, 56, 6, 155, 88, 26, 143, 121, 233, 112, 208, 194, 173, 168, 80, 117, 132, 72,
    202, 252, 218, 138, 214, 84, 66, 36, 191, 152, 149, 249, 227, 94, 181, 21, 104, 97, 40, 186,
    223, 76, 241, 47, 129, 230, 178, 63, 51, 238, 54, 16, 110, 24, 70, 166, 34, 136, 19, 247,
    45, 184, 14, 61, 245, 164, 57, 59, 7, 158, 156, 157, 89, 159, 27, 8, 144, 9, 122, 28,
    234, 160, 113, 90, 209, 29, 195, 123, 174, 10, 169, 145, 81, 91, 118, 114, 133, 161, 73, 235,
    203, 124, 253, 196, 219, 30, 139, 210, 215, 146, 85, 170, 67, 11, 37, 175, 192, 115, 153, 119,
    150, 92, 250, 82, 228, 236, 95, 74, 182, 162, 22, 134, 105, 197, 98, 254, 41, 125, 187, 204,
    224, 211, 77, 140, 242, 31, 48, 220, 130, 171, 231, 86, 179, 147, 64, 216, 52, 176, 239, 38,
    55, 12, 17, 68, 111, 120, 25, 154, 71, 116, 167, 193, 35, 83, 137, 251, 20, 93, 248, 151,
    46, 75, 185, 96, 15, 237, 62, 229, 246, 135, 165, 23, 58, 163, 60, 183 };

static const u8 GFC256_EXP_TABLE[512 * 2 + 1] = {
    1, 2, 4, 8, 16, 32, 64, 128, 135, 137, 149, 173, 221, 61, 122, 244, 111, 222, 59, 118,
    236, 95, 190, 251, 113, 226, 67, 134, 139, 145, 165, 205, 29, 58, 116, 232, 87, 174, 219, 49,
    98, 196, 15, 30, 60, 120, 240, 103, 206, 27, 54, 108, 216, 55, 110, 220, 63, 126, 252, 127,
    254, 123, 246, 107, 214, 43, 86, 172, 223, 57, 114, 228, 79, 158, 187, 241, 101, 202, 19, 38,
    76, 152, 183, 233, 85, 170, 211, 33, 66, 132, 143, 153, 181, 237, 93, 186, 243, 97, 194, 3,
    6, 12, 24, 48, 96, 192, 7, 14, 28, 56, 112, 224, 71, 142, 155, 177, 229, 77, 154, 179,
    225, 69, 138, 147, 161, 197, 13, 26, 52, 104, 208, 39, 78, 156, 191, 249, 117, 234, 83, 166,
    203, 17, 34, 68, 136, 151, 169, 213, 45, 90, 180, 239, 89, 178, 227, 65, 130, 131, 129, 133,
    141, 157, 189, 253, 125, 250, 115, 230, 75, 150, 171, 209, 37, 74, 148, 175, 217, 53, 106, 212,
    47, 94, 188, 255, 121, 242, 99, 198, 11, 22, 44, 88, 176, 231, 73, 146, 163, 193, 5, 10,
    20, 40, 80, 160, 199, 9, 18, 36, 72, 144, 167, 201, 21, 42, 84, 168, 215, 41, 82, 164,
    207, 25, 50, 100, 200, 23, 46, 92, 184, 247, 105, 210, 35, 70, 140, 159, 185, 245, 109, 218,
    51, 102, 204, 31, 62, 124, 248, 119, 238, 91, 182, 235, 81, 162, 195, 1, 2, 4, 8, 16,
    32, 64, 128, 135, 137, 149, 173, 221, 61, 122, 244, 111, 222, 59, 118, 236, 95, 190, 251, 113,
    226, 67, 134, 139, 145, 165, 205, 29, 58, 116, 232, 87, 174, 219, 49, 98, 196, 15, 30, 60,
    120, 240, 103, 206, 27, 54, 108, 216, 55, 110, 220, 63, 126, 252, 127, 254, 123, 246, 107, 214,
    43, 86, 172, 223, 57, 114, 228, 79, 158, 187, 241, 101, 202, 19, 38, 76, 152, 183, 233, 85,
    170, 211, 33, 66, 132, 143, 153, 181, 237, 93, 186, 243, 97, 194, 3, 6, 12, 24, 48, 96,
    192, 7, 14, 28, 56, 112, 224, 71, 142, 155, 177, 229, 77, 154, 179, 225, 69, 138, 147, 161,
    197, 13, 26, 52, 104, 208, 39, 78, 156, 191, 249, 117, 234, 83, 166, 203, 17, 34, 68, 136,
    151, 169, 213, 45, 90, 180, 239, 89, 178, 227, 65, 130, 131, 129, 133, 141, 157, 189, 253, 125,
    250, 115, 230, 75, 150, 171, 209, 37, 74, 148, 175, 217, 53, 106, 212, 47, 94, 188, 255, 121,
    242, 99, 198, 11, 22, 44, 88, 176, 231, 73, 146, 163, 193, 5, 10, 20, 40, 80, 160, 199,
    9, 18, 36, 72, 144, 167, 201, 21, 42, 84, 168, 215, 41, 82, 164, 207, 25, 50, 100, 200,
    23, 46, 92, 184, 247, 105, 210, 35, 70, 140, 159, 185, 245, 109, 218, 51, 102, 204, 31, 62,
    124, 248, 119, 238, 91, 182, 235, 81, 162, 195, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static const u8 GFC256_INV_TABLE[256] = {
    0, 1, 195, 130, 162, 126, 65, 90, 81, 54, 63, 172, 227, 104, 45, 42, 235, 155, 27, 53,
    220, 30, 86, 165, 178, 116, 52, 18, 213, 100, 21, 221, 182, 75, 142, 251, 206, 233, 217, 161,
    110, 219, 15, 44, 43, 14, 145, 241, 89, 215, 58, 244, 26, 19, 9, 80, 169, 99, 50, 245,
    201, 204, 173, 10, 91, 6, 230, 247, 71, 191, 190, 68, 103, 123, 183, 33, 175, 83, 147, 255,
    55, 8, 174, 77, 196, 209, 22, 164, 214, 48, 7, 64, 139, 157, 187, 140, 239, 129, 168, 57,
    29, 212, 122, 72, 13, 226, 202, 176, 199, 222, 40, 218, 151, 210, 242, 132, 25, 179, 185, 135,
    167, 228, 102, 73, 149, 153, 5, 163, 238, 97, 3, 194, 115, 243, 184, 119, 224, 248, 156, 92,
    95, 186, 34, 250, 240, 46, 254, 78, 152, 124, 211, 112, 148, 125, 234, 17, 138, 93, 188, 236,
    216, 39, 4, 127, 87, 23, 229, 120, 98, 56, 171, 170, 11, 62, 82, 76, 107, 203, 24, 117,
    192, 253, 32, 74, 134, 118, 141, 94, 158, 237, 70, 69, 180, 252, 131, 2, 84, 208, 223, 108,
    205, 60, 106, 177, 61, 200, 36, 232, 197, 85, 113, 150, 101, 28, 88, 49, 160, 38, 111, 41,
    20, 31, 109, 198, 136, 249, 105, 12, 121, 166, 66, 246, 207, 37, 154, 16, 159, 189, 128, 96,
    144, 47, 114, 133, 51, 59, 231, 67, 137, 225, 143, 35, 193, 181, 146, 79 };

u8 * CAT_RESTRICT GFC256_MUL_TABLE = 0;
u8 * CAT_RESTRICT GFC256_DIV_TABLE = 0;

static void GFC256Init()
{
    if (GFC256_MUL_TABLE)
    {
        return;
    }

    // Allocate table memory 65KB x 2
    GFC256_MUL_TABLE = new u8[256 * 256 * 2];
    GFC256_DIV_TABLE = GFC256_MUL_TABLE + 256 * 256;

    u8 *m = GFC256_MUL_TABLE, *d = GFC256_DIV_TABLE;

    // Unroll y = 0 subtable
    for (int x = 0; x < 256; ++x)
    {
        m[x] = d[x] = 0;
    }

    // For each other y value,
    for (int y = 1; y < 256; ++y)
    {
        // Calculate log(y) for mult and 255 - log(y) for div
        const u8 log_y = (u8)GFC256_LOG_TABLE[y];
        const u8 log_yn = 255 - log_y;

        // Next subtable
        m += 256;
        d += 256;

        // Unroll x = 0
        m[0] = 0;
        d[0] = 0;

        // Calculate x * y, x / y
        for (int x = 1; x < 256; ++x)
        {
            int log_x = GFC256_LOG_TABLE[x];

            m[x] = GFC256_EXP_TABLE[log_x + log_y];
            d[x] = GFC256_EXP_TABLE[log_x + log_yn];
        }
    }
}

// return x * y in GF(256)
// For repeated multiplication by a constant, it is faster to put the constant in y.
static CAT_INLINE u8 GFC256Multiply(u8 x, u8 y)
{
    return GFC256_MUL_TABLE[((u32)y << 8) + x];
}

// return x / y in GF(256)
// Memory-access optimized for constant divisors in y.
static CAT_INLINE u8 GFC256Divide(u8 x, u8 y)
{
    return GFC256_DIV_TABLE[((u32)y << 8) + x];
}


// a_ij = 1 / (x_i - y_j)

u8 calc_a_ij(u8 x_i, u8 y_j)
{
    return GFC256_INV_TABLE[x_i ^ y_j];
}

u8 matrix_get(u8* matrix, int width, int row, int col)
{
    return matrix[row * width + col];
}

void matrix_row_mult(u8* matrix, int width, int row, u8 x)
{
    for (int i = 0; i < width; ++i)
    {
        matrix[row * width + i] = GFC256Multiply(matrix[row * width + i], x);
    }
}

// dest += src * x
void matrix_row_muladd(u8* matrix, int width, int src_row, int dest_row, u8 x)
{
    for (int i = 0; i < width; ++i)
    {
        matrix[dest_row * width + i] ^= GFC256Multiply(matrix[src_row * width + i], x);
    }
}

void matrix_swap(u8* matrix, int width, int row_a, int row_b)
{
    for (int i = 0; i < width; ++i)
    {
        u8 x = matrix[width * row_a + i];
        u8 y = matrix[width * row_b + i];
        matrix[width * row_a + i] = y;
        matrix[width * row_b + i] = x;
    }
}

bool gaussian_elimination(u8* matrix, int width)
{
    for (int i = 0; i < width; ++i)
    {
        bool found = false;
        for (int j = i; j < width; ++j)
        {
            u8 sym = matrix_get(matrix, width, j, i);
            if (sym != 0)
            {
                // found pivot!
                matrix_swap(matrix, width, i, j);

                matrix_row_mult(matrix, width, i, GFC256_INV_TABLE[sym]);

                for (int k = 0; k < width; ++k)
                {
                    if (k != i)
                    {
                        matrix_row_muladd(matrix, width, i, k, matrix_get(matrix, width, k, i));
                    }
                }

                found = true;
                break;
            }
        }
        if (!found)
        {
            return false;
        }
    }
    return true;
}

void generate_cauchy(u8* matrix, int width)
{
    u8 X[128], Y[128];
    for (int i = 0; i < 128; ++i)
    {
        X[i] = (u8)i;
        Y[i] = 128 + (u8)i;
    }

    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            matrix[width * i + j] = calc_a_ij(X[i], Y[j]);
        }
    }
}

#include <stdlib.h>
#include <iostream>
#include <vector>
using namespace std;

class Decoder
{
    struct Data
    {
        int group;
        int id;
        int t;
    };

    vector<Data> DataHistory;

    struct Dupe
    {
        vector<int> unknowns;
        int id;
        int t;
    };

    vector<Dupe> DupeHistory;

    void ExpireHistory(int t_oldest)
    {
        bool modded;
        do
        {
            modded = false;
            for (int i = 0; i < (int)DataHistory.size(); ++i)
            {
                if (DataHistory[i].t < t_oldest)
                {
                    DataHistory.erase(DataHistory.begin() + i);
                    modded = true;
                    break;
                }
            }
            for (int i = 0; i < (int)DupeHistory.size(); ++i)
            {
                if (DupeHistory[i].t < t_oldest)
                {
                    DupeHistory.erase(DupeHistory.begin() + i);
                    modded = true;
                    break;
                }
            }
        } while (modded);
    }

    void eliminateUnknown(int sym_id)
    {
        for (int i = 0; i < (int)DupeHistory.size(); ++i)
        {
            Dupe& dupe = DupeHistory[i];

            for (int j = 0; j < (int)dupe.unknowns.size(); ++j)
            {
            }
        }
    }

    int latest_sym_id, current_group;

    void updateLastestSymId(int sym_id)
    {
        // Assumes no re-ordering

        if (latest_sym_id < sym_id)
        {
            // group increments when symbol id rolls over
            current_group++;
        }

        latest_sym_id = sym_id;
    }

public:
    Decoder()
    {
        current_group = 0;
        latest_sym_id = 0;
    }
    ~Decoder()
    {
    }

    void OnData(int sym_id, int t_sent)
    {
        updateLastestSymId(sym_id);

        Data d;
        d.group = current_group;
        d.id = sym_id;
        d.t = t_sent;

        history.push_back(d);
    }

    void OnErasure(int start_sym_id, int stop_sym_id, int dupe_id, int t_sent)
    {
    }

    void PrintMatrix()
    {
        int rows = 0, cols = 0;

        vector<int> cols;

        for (int i = 0; i < (int)DupeHistory.size(); ++i)
        {
            Dupe& dupe = DupeHistory[i];

            dupe.unknowns
        }
    }
};


struct Symbol
{
    int id;
    int t;
    bool lost;
};

struct Dupe
{
    int start_sym_id, stop_sym_id;
    int dupe_id;
    int t;
    bool lost;
    vector<int> lost_syms;
};

void traffic_sim()
{
    int history_ms = 250;
    int send_ms = 20;
    int pkt_max = 20;
    int dupe_ms = 20;
    int max_sym_id = 100;
    int max_dupe_id = 28;

    Decoder decoder;

    srand(0);

    vector<Symbol> history;
    vector<Dupe> dupes;

    int sym_id = 0;
    int dupe_id = 0;
    int t = 0;

    int last_dupe_ms = 0;

    for (;;)
    {
        int sent = rand() % pkt_max;

        int loss_prob = 3;

        for (int i = 0; i < sent; ++i)
        {
            bool lost = (rand() % 100) < loss_prob;
            if (lost)
            {
                if (loss_prob < 50)
                {
                    loss_prob = 50;
                }
                else
                {
                    loss_prob = 3;
                }
            }
            else
            {
                loss_prob = 3;
            }

            Symbol s;
            s.id = sym_id;
            s.t = t;
            s.lost = lost;

            // Erase any duplicate sym_ids from history
            bool modded;
            do
            {
                modded = false;
                for (int i = 0; i < (int)history.size(); ++i)
                {
                    if (history[i].id == sym_id)
                    {
                        history.erase(history.begin() + i);
                        modded = true;
                        break;
                    }
                }
            } while (modded);

            decoder.OnData(sym_id, t);

            history.push_back(s);

            sym_id++;
            if (sym_id >= max_sym_id)
            {
                sym_id = 0;
            }
        }

        // If simulation time indicates another dupe should be sent,
        if (t - last_dupe_ms >= dupe_ms)
        {
            // Erase everything older than the history window
            bool modded;
            do 
            {
                modded = false;
                for (int i = 0; i < (int)history.size(); ++i)
                {
                    if (history[i].t < t - history_ms)
                    {
                        history.erase(history.begin() + i);
                        modded = true;
                        break;
                    }
                }
            } while (modded);

            if (history.size() > 0)
            {
                //for (int dupe_i = 0; dupe_i < 2; ++dupe_i)
                {
                    bool lost = (rand() % 100) < loss_prob;
                    if (lost)
                    {
                        if (loss_prob < 50)
                        {
                            loss_prob = 50;
                        }
                        else
                        {
                            loss_prob = 3;
                        }
                    }
                    else
                    {
                        loss_prob = 3;
                    }

                    last_dupe_ms = t;

                    Dupe d;
                    d.dupe_id = dupe_id;
                    d.t = t;
                    for (int i = 0; i < (int)history.size(); ++i)
                    {
                        if (history[i].lost)
                        {
                            d.lost_syms.push_back(history[i].id);
                        }
                    }
                    d.start_sym_id = history[0].id;
                    d.stop_sym_id = history[history.size() - 1].id;
                    d.lost = lost;
                    dupes.push_back(d);

                    decoder.OnErasure(d.start_sym_id, d.stop_sym_id, dupe_id, t);

                    dupe_id++;
                    if (dupe_id >= max_dupe_id)
                    {
                        dupe_id = 0;
                    }

                    for (int i = 0; i < (int)dupes.size(); ++i)
                    {
                        if (dupes[i].lost_syms.size() > 0)
                        {
                            cout << dupes[i].dupe_id << " : ";
                            for (int j = 0; j < (int)dupes[i].lost_syms.size(); ++j)
                            {
                                cout << dupes[i].lost_syms[j] << " ";
                            }
                            cout << endl;
                        }
                    }

                    dupes.clear();
                }
            }
        }

        t += send_ms;
    }
}

int main()
{
    GFC256Init();

    traffic_sim();

    /*
          0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
        8 x x x x x
        9   x x x x x
        10    x x x x x
        11      x x x x x
        8         x x x x x
        9           x x x x x
        10            x x x x x
        11              x x x x x

        Upper triangular matrices are all invertible.

        For lower triangular matrices plus a band to the right of the diagonal that is T-wide at the top:
            x x 0 0
            x x x 0
            x x x x
        These are invertible if |X| >= T and |Y| >= T.

            This seems to imply that T restricts the history of packets.
    */
    u8 matrix[128 * 128];
    int width = 128;

    u8 X[128], Y[128];
    for (int i = 0; i < 128; ++i)
    {
        X[i] = (u8)i;
        Y[i] = 128 + (u8)i;
    }

    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            matrix[width * i + j] = calc_a_ij(X[i % 6], Y[j % 6]);
            if (j - 5 > i)
            {
                matrix[width * i + j] = 0;
            }
        }
    }

    bool result = gaussian_elimination(matrix, width);

    return 0;
}
