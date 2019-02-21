
##
def to_spglib(strc):
    """
    converting the pylada structure object
    to the format spglib can digest
    """

    import numpy as np
    from pylada import periodic_table

    cell      = np.transpose(strc.cell)*strc.scale

    positions = [atom.pos*strc.scale for atom in strc]
    # converting into fractional coordinates
    positions = [ np.dot(np.linalg.inv(np.transpose(cell)),pos) for pos in positions ]

    symbols   = [periodic_table.find(symbol=atom.type).atomic_number for atom in strc]

    return (cell, positions, symbols)


##
def from_spglib(strc):
    """
    converting the pylada structure object
    from the spglib format to pylada
    """

    import numpy as np
    from pylada import periodic_table
    from pylada.crystal import Structure

    out_s = Structure()
    out_s.scale = 1.

    cell      = strc[0]
    positions = strc[1]
    symbols   = strc[2]

    out_s.cell = np.transpose(cell)

    for ii in range(len(positions)):
        pp=np.dot(np.transpose(cell),positions[ii])
        ss=periodic_table.symbols[symbols[ii]-1]
        out_s.add_atom(pp[0],pp[1],pp[2],ss)

    return out_s
##
