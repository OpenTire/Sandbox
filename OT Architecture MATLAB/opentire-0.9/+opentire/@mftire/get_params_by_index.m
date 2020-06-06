%% get_params_by_index - Get parameters by index

% part of mftire 1.1.0
function values = get_params_by_index(t, indices)
  %assert(all(indices >= 1 & indices <= t.PARAM_INFO.param_count));

  values = zeros(size(indices));
  for i=1:length(indices)
    switch indices(i)
    case 1
        value = t.FITTYP;
    case 2
        value = t.PLUS;
    case 3
        value = t.TYRESIDE;
    case 4
        value = t.LONGVL;
    case 5
        value = t.VXLOW;
    case 6
        value = t.ROAD_INCREMENT;
    case 7
        value = t.ROAD_DIRECTION;
    case 8
        value = t.PROPERTY_FILE_FORMAT;
    case 9
        value = t.USE_MODE;
    case 10
        value = t.HMAX_LOCAL;
    case 11
        value = t.TIME_SWITCH_INTEG;
    case 12
        value = t.USER_SUB_ID;
    case 13
        value = t.N_TIRE_STATES;
    case 14
        value = t.UNLOADED_RADIUS;
    case 15
        value = t.WIDTH;
    case 16
        value = t.RIM_RADIUS;
    case 17
        value = t.RIM_WIDTH;
    case 18
        value = t.ASPECT_RATIO;
    case 19
        value = t.INFLPRES;
    case 20
        value = t.NOMPRES;
    case 21
        value = t.MASS;
    case 22
        value = t.IXX;
    case 23
        value = t.IYY;
    case 24
        value = t.BELT_MASS;
    case 25
        value = t.BELT_IXX;
    case 26
        value = t.BELT_IYY;
    case 27
        value = t.GRAVITY;
    case 28
        value = t.FNOMIN;
    case 29
        value = t.VERTICAL_STIFFNESS;
    case 30
        value = t.VERTICAL_DAMPING;
    case 31
        value = t.MC_CONTOUR_A;
    case 32
        value = t.MC_CONTOUR_B;
    case 33
        value = t.BREFF;
    case 34
        value = t.DREFF;
    case 35
        value = t.FREFF;
    case 36
        value = t.Q_RE0;
    case 37
        value = t.Q_V1;
    case 38
        value = t.Q_V2;
    case 39
        value = t.Q_FZ2;
    case 40
        value = t.Q_FCX;
    case 41
        value = t.Q_FCY;
    case 42
        value = t.Q_FCY2;
    case 43
        value = t.Q_CAM;
    case 44
        value = t.Q_CAM1;
    case 45
        value = t.Q_CAM2;
    case 46
        value = t.Q_CAM3;
    case 47
        value = t.Q_FYS1;
    case 48
        value = t.Q_FYS2;
    case 49
        value = t.Q_FYS3;
    case 50
        value = t.Q_FZ3;
    case 51
        value = t.PFZ1;
    case 52
        value = t.BOTTOM_OFFST;
    case 53
        value = t.BOTTOM_STIFF;
    case 54
        value = t.LONGITUDINAL_STIFFNESS;
    case 55
        value = t.LATERAL_STIFFNESS;
    case 56
        value = t.YAW_STIFFNESS;
    case 57
        value = t.FREQ_LONG;
    case 58
        value = t.FREQ_LAT;
    case 59
        value = t.FREQ_YAW;
    case 60
        value = t.FREQ_WINDUP;
    case 61
        value = t.DAMP_LONG;
    case 62
        value = t.DAMP_LAT;
    case 63
        value = t.DAMP_YAW;
    case 64
        value = t.DAMP_WINDUP;
    case 65
        value = t.DAMP_RESIDUAL;
    case 66
        value = t.DAMP_VLOW;
    case 67
        value = t.Q_BVX;
    case 68
        value = t.Q_BVT;
    case 69
        value = t.PCFX1;
    case 70
        value = t.PCFX2;
    case 71
        value = t.PCFX3;
    case 72
        value = t.PCFY1;
    case 73
        value = t.PCFY2;
    case 74
        value = t.PCFY3;
    case 75
        value = t.PCMZ1;
    case 76
        value = t.Q_RA1;
    case 77
        value = t.Q_RA2;
    case 78
        value = t.Q_RB1;
    case 79
        value = t.Q_RB2;
    case 80
        value = t.ELLIPS_SHIFT;
    case 81
        value = t.ELLIPS_LENGTH;
    case 82
        value = t.ELLIPS_HEIGHT;
    case 83
        value = t.ELLIPS_ORDER;
    case 84
        value = t.ELLIPS_MAX_STEP;
    case 85
        value = t.ELLIPS_NWIDTH;
    case 86
        value = t.ELLIPS_NLENGTH;
    case 87
        value = t.ENV_C1;
    case 88
        value = t.ENV_C2;
    case 89
        value = t.Q_A1;
    case 90
        value = t.Q_A2;
    case 91
        value = t.PRESMIN;
    case 92
        value = t.PRESMAX;
    case 93
        value = t.FZMIN;
    case 94
        value = t.FZMAX;
    case 95
        value = t.KPUMIN;
    case 96
        value = t.KPUMAX;
    case 97
        value = t.ALPMIN;
    case 98
        value = t.ALPMAX;
    case 99
        value = t.CAMMIN;
    case 100
        value = t.CAMMAX;
    case 101
        value = t.LFZO;
    case 102
        value = t.LCX;
    case 103
        value = t.LMUX;
    case 104
        value = t.LEX;
    case 105
        value = t.LKX;
    case 106
        value = t.LHX;
    case 107
        value = t.LVX;
    case 108
        value = t.LCY;
    case 109
        value = t.LMUY;
    case 110
        value = t.LEY;
    case 111
        value = t.LKY;
    case 112
        value = t.LKYC;
    case 113
        value = t.LKZC;
    case 114
        value = t.LHY;
    case 115
        value = t.LVY;
    case 116
        value = t.LTR;
    case 117
        value = t.LRES;
    case 118
        value = t.LXAL;
    case 119
        value = t.LYKA;
    case 120
        value = t.LVYKA;
    case 121
        value = t.LS;
    case 122
        value = t.LMX;
    case 123
        value = t.LVMX;
    case 124
        value = t.LMY;
    case 125
        value = t.LMP;
    case 126
        value = t.LSGKP;
    case 127
        value = t.LSGAL;
    case 128
        value = t.LMUV;
    case 129
        value = t.LAMU;
    case 130
        value = t.PCX1;
    case 131
        value = t.PDX1;
    case 132
        value = t.PDX2;
    case 133
        value = t.PDX3;
    case 134
        value = t.PEX1;
    case 135
        value = t.PEX2;
    case 136
        value = t.PEX3;
    case 137
        value = t.PEX4;
    case 138
        value = t.PKX1;
    case 139
        value = t.PKX2;
    case 140
        value = t.PKX3;
    case 141
        value = t.PHX1;
    case 142
        value = t.PHX2;
    case 143
        value = t.PVX1;
    case 144
        value = t.PVX2;
    case 145
        value = t.RBX1;
    case 146
        value = t.RBX2;
    case 147
        value = t.RBX3;
    case 148
        value = t.RCX1;
    case 149
        value = t.REX1;
    case 150
        value = t.REX2;
    case 151
        value = t.RHX1;
    case 152
        value = t.PPX1;
    case 153
        value = t.PPX2;
    case 154
        value = t.PPX3;
    case 155
        value = t.PPX4;
    case 156
        value = t.PTX1;
    case 157
        value = t.PTX2;
    case 158
        value = t.PTX3;
    case 159
        value = t.QSX1;
    case 160
        value = t.QSX2;
    case 161
        value = t.QSX3;
    case 162
        value = t.QSX4;
    case 163
        value = t.QSX5;
    case 164
        value = t.QSX6;
    case 165
        value = t.QSX7;
    case 166
        value = t.QSX8;
    case 167
        value = t.QSX9;
    case 168
        value = t.QSX10;
    case 169
        value = t.QSX11;
    case 170
        value = t.QSX12;
    case 171
        value = t.QSX13;
    case 172
        value = t.QSX14;
    case 173
        value = t.PPMX1;
    case 174
        value = t.PCY1;
    case 175
        value = t.PDY1;
    case 176
        value = t.PDY2;
    case 177
        value = t.PDY3;
    case 178
        value = t.PEY1;
    case 179
        value = t.PEY2;
    case 180
        value = t.PEY3;
    case 181
        value = t.PEY4;
    case 182
        value = t.PEY5;
    case 183
        value = t.PKY1;
    case 184
        value = t.PKY2;
    case 185
        value = t.PKY3;
    case 186
        value = t.PKY4;
    case 187
        value = t.PKY5;
    case 188
        value = t.PKY6;
    case 189
        value = t.PKY7;
    case 190
        value = t.PHY1;
    case 191
        value = t.PHY2;
    case 192
        value = t.PHY3;
    case 193
        value = t.PVY1;
    case 194
        value = t.PVY2;
    case 195
        value = t.PVY3;
    case 196
        value = t.PVY4;
    case 197
        value = t.RBY1;
    case 198
        value = t.RBY2;
    case 199
        value = t.RBY3;
    case 200
        value = t.RBY4;
    case 201
        value = t.RCY1;
    case 202
        value = t.REY1;
    case 203
        value = t.REY2;
    case 204
        value = t.RHY1;
    case 205
        value = t.RHY2;
    case 206
        value = t.RVY1;
    case 207
        value = t.RVY2;
    case 208
        value = t.RVY3;
    case 209
        value = t.RVY4;
    case 210
        value = t.RVY5;
    case 211
        value = t.RVY6;
    case 212
        value = t.PPY1;
    case 213
        value = t.PPY2;
    case 214
        value = t.PPY3;
    case 215
        value = t.PPY4;
    case 216
        value = t.PPY5;
    case 217
        value = t.PTY1;
    case 218
        value = t.PTY2;
    case 219
        value = t.QSY1;
    case 220
        value = t.QSY2;
    case 221
        value = t.QSY3;
    case 222
        value = t.QSY4;
    case 223
        value = t.QSY5;
    case 224
        value = t.QSY6;
    case 225
        value = t.QSY7;
    case 226
        value = t.QSY8;
    case 227
        value = t.QBZ1;
    case 228
        value = t.QBZ2;
    case 229
        value = t.QBZ3;
    case 230
        value = t.QBZ4;
    case 231
        value = t.QBZ5;
    case 232
        value = t.QBZ9;
    case 233
        value = t.QBZ10;
    case 234
        value = t.QCZ1;
    case 235
        value = t.QDZ1;
    case 236
        value = t.QDZ2;
    case 237
        value = t.QDZ3;
    case 238
        value = t.QDZ4;
    case 239
        value = t.QDZ6;
    case 240
        value = t.QDZ7;
    case 241
        value = t.QDZ8;
    case 242
        value = t.QDZ9;
    case 243
        value = t.QDZ10;
    case 244
        value = t.QDZ11;
    case 245
        value = t.QEZ1;
    case 246
        value = t.QEZ2;
    case 247
        value = t.QEZ3;
    case 248
        value = t.QEZ4;
    case 249
        value = t.QEZ5;
    case 250
        value = t.QHZ1;
    case 251
        value = t.QHZ2;
    case 252
        value = t.QHZ3;
    case 253
        value = t.QHZ4;
    case 254
        value = t.SSZ1;
    case 255
        value = t.SSZ2;
    case 256
        value = t.SSZ3;
    case 257
        value = t.SSZ4;
    case 258
        value = t.PPZ1;
    case 259
        value = t.PPZ2;
    case 260
        value = t.PDXP1;
    case 261
        value = t.PDXP2;
    case 262
        value = t.PDXP3;
    case 263
        value = t.PKYP1;
    case 264
        value = t.PDYP1;
    case 265
        value = t.PDYP2;
    case 266
        value = t.PDYP3;
    case 267
        value = t.PDYP4;
    case 268
        value = t.PHYP1;
    case 269
        value = t.PHYP2;
    case 270
        value = t.PHYP3;
    case 271
        value = t.PHYP4;
    case 272
        value = t.PECP1;
    case 273
        value = t.PECP2;
    case 274
        value = t.QDTP1;
    case 275
        value = t.QCRP1;
    case 276
        value = t.QCRP2;
    case 277
        value = t.QBRP1;
    case 278
        value = t.QDRP1;
    case 279
        value = t.QDRP2;
    otherwise
    end
    values(i) = value;
  end
end
