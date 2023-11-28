def calc_forces(filenum = filenumber):
    pos = positions[filenum]
    dx, dy, dz, rik = relative_pos(pos)
    f = np.zeros((n,3))
    nghbrs_count, nghbrs_indices = neighbors_list()
    for k in range(n):
        for j in range(int(nghbrs_count[k])):
            i = nghbrs_indices[k][j]
            drik = rik[i][k]
            if drik != 0:
                c1 = 24*epsilon*(sigma**6/drik**8)*((2*sigma**6)/drik**6 - 1)
                f[k,0] += c1*dx[i][k]
                f[k,1] += c1*dy[i][k]
                f[k,2] += c1*dz[i][k]
    return f