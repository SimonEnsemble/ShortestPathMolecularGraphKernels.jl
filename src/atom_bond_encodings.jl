# for speed, we encode the bond and atom types at UInt8's.

"""
    BondType(bond_order::Int)

an efficient `UInt8` encoding of a bond type.

# example
```julia
bt = BondType(0) # BOND_AROMATIC
bt = BondType(1) # BOND_SINGLE
bt = BondType(2) # BOND_DOUBLE
bt = BondType(3) # BOND_TRIPLE
```
"""
@enum BondType::UInt8 begin
    BOND_AROMATIC = 0
    BOND_SINGLE = 1
    BOND_DOUBLE = 2
    BOND_TRIPLE = 3
end

# generate via PeriodicTable.jl:
# for element in elements
#   println("ATOM_" * element.symbol * " = " * "$(element.number)")
# end
"""
    AtomType(atomic_number::Int)

an efficient `UInt8` encoding of an atom type. 

# example
```julia
at = AtomType(6) # ATOM_C
```
"""
@enum AtomType::UInt8 begin
    ATOM_H = 1
    ATOM_He = 2
    ATOM_Li = 3
    ATOM_Be = 4
    ATOM_B = 5
    ATOM_C = 6
    ATOM_N = 7
    ATOM_O = 8
    ATOM_F = 9
    ATOM_Ne = 10
    ATOM_Na = 11
    ATOM_Mg = 12
    ATOM_Al = 13
    ATOM_Si = 14
    ATOM_P = 15
    ATOM_S = 16
    ATOM_Cl = 17
    ATOM_Ar = 18
    ATOM_K = 19
    ATOM_Ca = 20
    ATOM_Sc = 21
    ATOM_Ti = 22
    ATOM_V = 23
    ATOM_Cr = 24
    ATOM_Mn = 25
    ATOM_Fe = 26
    ATOM_Co = 27
    ATOM_Ni = 28
    ATOM_Cu = 29
    ATOM_Zn = 30
    ATOM_Ga = 31
    ATOM_Ge = 32
    ATOM_As = 33
    ATOM_Se = 34
    ATOM_Br = 35
    ATOM_Kr = 36
    ATOM_Rb = 37
    ATOM_Sr = 38
    ATOM_Y = 39
    ATOM_Zr = 40
    ATOM_Nb = 41
    ATOM_Mo = 42
    ATOM_Tc = 43
    ATOM_Ru = 44
    ATOM_Rh = 45
    ATOM_Pd = 46
    ATOM_Ag = 47
    ATOM_Cd = 48
    ATOM_In = 49
    ATOM_Sn = 50
    ATOM_Sb = 51
    ATOM_Te = 52
    ATOM_I = 53
    ATOM_Xe = 54
    ATOM_Cs = 55
    ATOM_Ba = 56
    ATOM_La = 57
    ATOM_Ce = 58
    ATOM_Pr = 59
    ATOM_Nd = 60
    ATOM_Pm = 61
    ATOM_Sm = 62
    ATOM_Eu = 63
    ATOM_Gd = 64
    ATOM_Tb = 65
    ATOM_Dy = 66
    ATOM_Ho = 67
    ATOM_Er = 68
    ATOM_Tm = 69
    ATOM_Yb = 70
    ATOM_Lu = 71
    ATOM_Hf = 72
    ATOM_Ta = 73
    ATOM_W = 74
    ATOM_Re = 75
    ATOM_Os = 76
    ATOM_Ir = 77
    ATOM_Pt = 78
    ATOM_Au = 79
    ATOM_Hg = 80
    ATOM_Tl = 81
    ATOM_Pb = 82
    ATOM_Bi = 83
    ATOM_Po = 84
    ATOM_At = 85
    ATOM_Rn = 86
    ATOM_Fr = 87
    ATOM_Ra = 88
    ATOM_Ac = 89
    ATOM_Th = 90
    ATOM_Pa = 91
    ATOM_U = 92
    ATOM_Np = 93
    ATOM_Pu = 94
    ATOM_Am = 95
    ATOM_Cm = 96
    ATOM_Bk = 97
    ATOM_Cf = 98
    ATOM_Es = 99
    ATOM_Fm = 100
    ATOM_Md = 101
    ATOM_No = 102
    ATOM_Lr = 103
    ATOM_Rf = 104
    ATOM_Db = 105
    ATOM_Sg = 106
    ATOM_Bh = 107
    ATOM_Hs = 108
    ATOM_Mt = 109
    ATOM_Ds = 110
    ATOM_Rg = 111
    ATOM_Cn = 112
    ATOM_Nh = 113
    ATOM_Fl = 114
    ATOM_Mc = 115
    ATOM_Lv = 116
    ATOM_Ts = 117
    ATOM_Og = 118
    ATOM_Uue = 119
end
