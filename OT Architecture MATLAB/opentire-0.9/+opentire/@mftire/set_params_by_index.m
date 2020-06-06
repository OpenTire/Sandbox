%% set_params_by_index - Set parameters by index

% part of mftire 1.1.0
function set_params_by_index(t, indices, values)
  %assert(all(indices >= 1 & indices <= t.PARAM_INFO.param_count));

  for i=1:length(indices)
    value = values(i);
    switch indices(i)
    case 1
        t.FITTYP = value;
    case 2
        t.PLUS = value;
    case 3
        t.TYRESIDE = value;
    case 4
        t.LONGVL = value;
    case 5
        t.VXLOW = value;
    case 6
        t.ROAD_INCREMENT = value;
    case 7
        t.ROAD_DIRECTION = value;
    case 8
        t.PROPERTY_FILE_FORMAT = value;
    case 9
        t.USE_MODE = value;
    case 10
        t.HMAX_LOCAL = value;
    case 11
        t.TIME_SWITCH_INTEG = value;
    case 12
        t.USER_SUB_ID = value;
    case 13
        t.N_TIRE_STATES = value;
    case 14
        t.UNLOADED_RADIUS = value;
    case 15
        t.WIDTH = value;
    case 16
        t.RIM_RADIUS = value;
    case 17
        t.RIM_WIDTH = value;
    case 18
        t.ASPECT_RATIO = value;
    case 19
        t.INFLPRES = value;
    case 20
        t.NOMPRES = value;
    case 21
        t.MASS = value;
    case 22
        t.IXX = value;
    case 23
        t.IYY = value;
    case 24
        t.BELT_MASS = value;
    case 25
        t.BELT_IXX = value;
    case 26
        t.BELT_IYY = value;
    case 27
        t.GRAVITY = value;
    case 28
        t.FNOMIN = value;
    case 29
        t.VERTICAL_STIFFNESS = value;
    case 30
        t.VERTICAL_DAMPING = value;
    case 31
        t.MC_CONTOUR_A = value;
    case 32
        t.MC_CONTOUR_B = value;
    case 33
        t.BREFF = value;
    case 34
        t.DREFF = value;
    case 35
        t.FREFF = value;
    case 36
        t.Q_RE0 = value;
    case 37
        t.Q_V1 = value;
    case 38
        t.Q_V2 = value;
    case 39
        t.Q_FZ2 = value;
    case 40
        t.Q_FCX = value;
    case 41
        t.Q_FCY = value;
    case 42
        t.Q_FCY2 = value;
    case 43
        t.Q_CAM = value;
    case 44
        t.Q_CAM1 = value;
    case 45
        t.Q_CAM2 = value;
    case 46
        t.Q_CAM3 = value;
    case 47
        t.Q_FYS1 = value;
    case 48
        t.Q_FYS2 = value;
    case 49
        t.Q_FYS3 = value;
    case 50
        t.Q_FZ3 = value;
    case 51
        t.PFZ1 = value;
    case 52
        t.BOTTOM_OFFST = value;
    case 53
        t.BOTTOM_STIFF = value;
    case 54
        t.LONGITUDINAL_STIFFNESS = value;
    case 55
        t.LATERAL_STIFFNESS = value;
    case 56
        t.YAW_STIFFNESS = value;
    case 57
        t.FREQ_LONG = value;
    case 58
        t.FREQ_LAT = value;
    case 59
        t.FREQ_YAW = value;
    case 60
        t.FREQ_WINDUP = value;
    case 61
        t.DAMP_LONG = value;
    case 62
        t.DAMP_LAT = value;
    case 63
        t.DAMP_YAW = value;
    case 64
        t.DAMP_WINDUP = value;
    case 65
        t.DAMP_RESIDUAL = value;
    case 66
        t.DAMP_VLOW = value;
    case 67
        t.Q_BVX = value;
    case 68
        t.Q_BVT = value;
    case 69
        t.PCFX1 = value;
    case 70
        t.PCFX2 = value;
    case 71
        t.PCFX3 = value;
    case 72
        t.PCFY1 = value;
    case 73
        t.PCFY2 = value;
    case 74
        t.PCFY3 = value;
    case 75
        t.PCMZ1 = value;
    case 76
        t.Q_RA1 = value;
    case 77
        t.Q_RA2 = value;
    case 78
        t.Q_RB1 = value;
    case 79
        t.Q_RB2 = value;
    case 80
        t.ELLIPS_SHIFT = value;
    case 81
        t.ELLIPS_LENGTH = value;
    case 82
        t.ELLIPS_HEIGHT = value;
    case 83
        t.ELLIPS_ORDER = value;
    case 84
        t.ELLIPS_MAX_STEP = value;
    case 85
        t.ELLIPS_NWIDTH = value;
    case 86
        t.ELLIPS_NLENGTH = value;
    case 87
        t.ENV_C1 = value;
    case 88
        t.ENV_C2 = value;
    case 89
        t.Q_A1 = value;
    case 90
        t.Q_A2 = value;
    case 91
        t.PRESMIN = value;
    case 92
        t.PRESMAX = value;
    case 93
        t.FZMIN = value;
    case 94
        t.FZMAX = value;
    case 95
        t.KPUMIN = value;
    case 96
        t.KPUMAX = value;
    case 97
        t.ALPMIN = value;
    case 98
        t.ALPMAX = value;
    case 99
        t.CAMMIN = value;
    case 100
        t.CAMMAX = value;
    case 101
        t.LFZO = value;
    case 102
        t.LCX = value;
    case 103
        t.LMUX = value;
    case 104
        t.LEX = value;
    case 105
        t.LKX = value;
    case 106
        t.LHX = value;
    case 107
        t.LVX = value;
    case 108
        t.LCY = value;
    case 109
        t.LMUY = value;
    case 110
        t.LEY = value;
    case 111
        t.LKY = value;
    case 112
        t.LKYC = value;
    case 113
        t.LKZC = value;
    case 114
        t.LHY = value;
    case 115
        t.LVY = value;
    case 116
        t.LTR = value;
    case 117
        t.LRES = value;
    case 118
        t.LXAL = value;
    case 119
        t.LYKA = value;
    case 120
        t.LVYKA = value;
    case 121
        t.LS = value;
    case 122
        t.LMX = value;
    case 123
        t.LVMX = value;
    case 124
        t.LMY = value;
    case 125
        t.LMP = value;
    case 126
        t.LSGKP = value;
    case 127
        t.LSGAL = value;
    case 128
        t.LMUV = value;
    case 129
        t.LAMU = value;
    case 130
        t.PCX1 = value;
    case 131
        t.PDX1 = value;
    case 132
        t.PDX2 = value;
    case 133
        t.PDX3 = value;
    case 134
        t.PEX1 = value;
    case 135
        t.PEX2 = value;
    case 136
        t.PEX3 = value;
    case 137
        t.PEX4 = value;
    case 138
        t.PKX1 = value;
    case 139
        t.PKX2 = value;
    case 140
        t.PKX3 = value;
    case 141
        t.PHX1 = value;
    case 142
        t.PHX2 = value;
    case 143
        t.PVX1 = value;
    case 144
        t.PVX2 = value;
    case 145
        t.RBX1 = value;
    case 146
        t.RBX2 = value;
    case 147
        t.RBX3 = value;
    case 148
        t.RCX1 = value;
    case 149
        t.REX1 = value;
    case 150
        t.REX2 = value;
    case 151
        t.RHX1 = value;
    case 152
        t.PPX1 = value;
    case 153
        t.PPX2 = value;
    case 154
        t.PPX3 = value;
    case 155
        t.PPX4 = value;
    case 156
        t.PTX1 = value;
    case 157
        t.PTX2 = value;
    case 158
        t.PTX3 = value;
    case 159
        t.QSX1 = value;
    case 160
        t.QSX2 = value;
    case 161
        t.QSX3 = value;
    case 162
        t.QSX4 = value;
    case 163
        t.QSX5 = value;
    case 164
        t.QSX6 = value;
    case 165
        t.QSX7 = value;
    case 166
        t.QSX8 = value;
    case 167
        t.QSX9 = value;
    case 168
        t.QSX10 = value;
    case 169
        t.QSX11 = value;
    case 170
        t.QSX12 = value;
    case 171
        t.QSX13 = value;
    case 172
        t.QSX14 = value;
    case 173
        t.PPMX1 = value;
    case 174
        t.PCY1 = value;
    case 175
        t.PDY1 = value;
    case 176
        t.PDY2 = value;
    case 177
        t.PDY3 = value;
    case 178
        t.PEY1 = value;
    case 179
        t.PEY2 = value;
    case 180
        t.PEY3 = value;
    case 181
        t.PEY4 = value;
    case 182
        t.PEY5 = value;
    case 183
        t.PKY1 = value;
    case 184
        t.PKY2 = value;
    case 185
        t.PKY3 = value;
    case 186
        t.PKY4 = value;
    case 187
        t.PKY5 = value;
    case 188
        t.PKY6 = value;
    case 189
        t.PKY7 = value;
    case 190
        t.PHY1 = value;
    case 191
        t.PHY2 = value;
    case 192
        t.PHY3 = value;
    case 193
        t.PVY1 = value;
    case 194
        t.PVY2 = value;
    case 195
        t.PVY3 = value;
    case 196
        t.PVY4 = value;
    case 197
        t.RBY1 = value;
    case 198
        t.RBY2 = value;
    case 199
        t.RBY3 = value;
    case 200
        t.RBY4 = value;
    case 201
        t.RCY1 = value;
    case 202
        t.REY1 = value;
    case 203
        t.REY2 = value;
    case 204
        t.RHY1 = value;
    case 205
        t.RHY2 = value;
    case 206
        t.RVY1 = value;
    case 207
        t.RVY2 = value;
    case 208
        t.RVY3 = value;
    case 209
        t.RVY4 = value;
    case 210
        t.RVY5 = value;
    case 211
        t.RVY6 = value;
    case 212
        t.PPY1 = value;
    case 213
        t.PPY2 = value;
    case 214
        t.PPY3 = value;
    case 215
        t.PPY4 = value;
    case 216
        t.PPY5 = value;
    case 217
        t.PTY1 = value;
    case 218
        t.PTY2 = value;
    case 219
        t.QSY1 = value;
    case 220
        t.QSY2 = value;
    case 221
        t.QSY3 = value;
    case 222
        t.QSY4 = value;
    case 223
        t.QSY5 = value;
    case 224
        t.QSY6 = value;
    case 225
        t.QSY7 = value;
    case 226
        t.QSY8 = value;
    case 227
        t.QBZ1 = value;
    case 228
        t.QBZ2 = value;
    case 229
        t.QBZ3 = value;
    case 230
        t.QBZ4 = value;
    case 231
        t.QBZ5 = value;
    case 232
        t.QBZ9 = value;
    case 233
        t.QBZ10 = value;
    case 234
        t.QCZ1 = value;
    case 235
        t.QDZ1 = value;
    case 236
        t.QDZ2 = value;
    case 237
        t.QDZ3 = value;
    case 238
        t.QDZ4 = value;
    case 239
        t.QDZ6 = value;
    case 240
        t.QDZ7 = value;
    case 241
        t.QDZ8 = value;
    case 242
        t.QDZ9 = value;
    case 243
        t.QDZ10 = value;
    case 244
        t.QDZ11 = value;
    case 245
        t.QEZ1 = value;
    case 246
        t.QEZ2 = value;
    case 247
        t.QEZ3 = value;
    case 248
        t.QEZ4 = value;
    case 249
        t.QEZ5 = value;
    case 250
        t.QHZ1 = value;
    case 251
        t.QHZ2 = value;
    case 252
        t.QHZ3 = value;
    case 253
        t.QHZ4 = value;
    case 254
        t.SSZ1 = value;
    case 255
        t.SSZ2 = value;
    case 256
        t.SSZ3 = value;
    case 257
        t.SSZ4 = value;
    case 258
        t.PPZ1 = value;
    case 259
        t.PPZ2 = value;
    case 260
        t.PDXP1 = value;
    case 261
        t.PDXP2 = value;
    case 262
        t.PDXP3 = value;
    case 263
        t.PKYP1 = value;
    case 264
        t.PDYP1 = value;
    case 265
        t.PDYP2 = value;
    case 266
        t.PDYP3 = value;
    case 267
        t.PDYP4 = value;
    case 268
        t.PHYP1 = value;
    case 269
        t.PHYP2 = value;
    case 270
        t.PHYP3 = value;
    case 271
        t.PHYP4 = value;
    case 272
        t.PECP1 = value;
    case 273
        t.PECP2 = value;
    case 274
        t.QDTP1 = value;
    case 275
        t.QCRP1 = value;
    case 276
        t.QCRP2 = value;
    case 277
        t.QBRP1 = value;
    case 278
        t.QDRP1 = value;
    case 279
        t.QDRP2 = value;
    otherwise
    end
  end
end
