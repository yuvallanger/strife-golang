// Auxilliary function to get me a python-like modulus.
// Used to get coordinates in a toroid board.

long SIGNAL = {signal};
long RECEPTOR = {receptor};
long COOPERATION = {cooperation};
long S_rad = {S_rad};
long C_rad = {C_rad};
long Nboard_row = {Nboard_row};
long Nboard_col = {Nboard_col};

long python_modulus(long index, long size)
{
    return ((index % size) + size);
}

// s_count() will only count for signallers.
// TODO - Add the "signal_kind" argument for Avigdor's model.
long s_count(long row_center, long col_center)
{
    long nh_row_i, nh_col_i;
    long count = 0;

    for (nh_row_i = -S_rad; nh_row_i <= S_rad; nh_row_i++)
    {
        for (nh_col_i = -S_rad; nh_col_i <= S_rad; nh_col_i++)
        {
            if (BOARD3(python_modulus(row_center + S_rad, Nboard_row),
                       python_modulus(col_center + S_rad, Nboard_col),
                       SIGNAL) == 1)
            {
                count++;
            }
        }
    }
    return count;
}

long is_c(long row_center, long col_center)
{
    long nh_row_i, nh_col_i;

    // A cell will make public goods if
    if (BOARD3(python_modulus(row_center, Nboard_row),
               python_modulus(col_center, Nboard_col),
               COOPERATION))
    {
        if (!BOARD3(python_modulus(row_center, Nboard_row), // The cell is a cooperator and
                    python_modulus(col_center, Nboard_col), //   does not have a receptor.
                    RECEPTOR) ||
            // The cell is a cooperator, has a working receptor and there is enough signal.
            (s_count(row_center, col_center) >= S_th))
        {
            return 1;
        }
    }
    else
    {
        return 0;
    }
}

long c_count(long row_center, long col_center)
{
    long nh_row_i, nh_col_i;
    long count = 0;

    for (nh_row_i = -C_rad; nh_row_i <= C_rad; nh_row_i++)
    {
        for (nh_col_i = -C_rad; nh_col_i <= C_rad; nh_col_i++)
        {
            if (is_c(row_center + nh_col_i,
                     col_center + nh_col_i))
            {
                count++;
            }
        }
    }

    return count;
}

long metabolism(long row, long col)
{
    return c_count(row, col) * (1 - {benefit}) * ({R_cost} * BOARD3(row, col, RECEPTOR) +
                                                  {S_cost} * BOARD3(row, col, SIGNAL) +
                                                  {C_cost} * BOARD3(row, col, COOPERATION));
}
