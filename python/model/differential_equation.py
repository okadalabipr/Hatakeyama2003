def diffeq(y,t,*args):

  v = [0]*35 #Rate equations

  v[1] = x[kf1]*y[R]*y[HRG] - x[kr1]*y[R_HRG]
  v[2] = x[kf2]*y[R_HRG]**2 - x[kr2]*y[R_HRG2]
  v[3] = x[kf3]*y[R_HRG2] - x[kr3]*y[RP]
  v[4] = x[V4]*y[RP]/(x[K4]+y[RP])
  v[5] = x[kf5]*y[RP]*y[Shc] - x[kr5]*y[R_Shc]
  v[6] = x[kf6]*y[R_Shc] - x[kr6]*y[R_ShP]
  v[7] = x[kf7]*y[R_ShP]*y[GS] - x[kr7]*y[R_ShGS]
  v[8] = x[kf8]*y[R_ShGS] - x[kr8]*y[ShGS]*y[RP]
  v[9] = x[kf9]*y[ShGS] - x[kr9]*y[GS]*y[ShP]
  v[10] = x[V10]*y[ShP]/(x[K10]+y[ShP])
  v[11] = x[kf11]*y[ShGS]*y[RasGDP]/(x[K11]+y[RasGDP])
  v[12] = x[V12]*y[RasGTP]/(x[K12]+y[RasGTP])
  v[13] = x[kf13]*y[RasGTP]*y[Raf]/(x[K13]+y[Raf])
  v[14] = x[kf14]*(y[Akt_PI_PP]+y[E])*y[Raf_act]/(x[K14]+y[Raf_act])
  v[15] = x[kf15]*y[Raf_act]*y[MEK]/(x[K15]*(1+y[MEKP]/x[K17])+y[MEK])
  v[16] = x[kf16]*y[PP2A]*y[MEKP]/(x[K16]*(1+y[MEKPP]/x[K18]+y[Akt_PI_P]/x[K31]+y[Akt_PI_PP]/x[K33])+y[MEKP])
  v[17] = x[kf17]*y[Raf_act]*y[MEKP]/(x[K17]*(1+y[MEK]/x[K15])+y[MEKP])
  v[18] = x[kf18]*y[PP2A]*y[MEKPP]/(x[K18]*(1+y[MEKP]/x[K16]+y[Akt_PI_P]/x[K31]+y[Akt_PI_PP]/x[K33])+y[MEKPP])
  v[19] = x[kf19]*y[MEKPP]*y[ERK]/(x[K19]*(1+y[ERKP]/x[K21])+y[ERK])
  v[20] = x[kf20]*y[MKP3]*y[ERKP]/(x[K20]*(1+y[ERKPP]/x[K22])+y[ERKP])
  v[21] = x[kf21]*y[MEKPP]*y[ERKP]/(x[K21]*(1+y[ERK]/x[K19])+y[ERKP])
  v[22] = x[kf22]*y[MKP3]*y[ERKPP]/(x[K22]*(1+y[ERKP]/x[K20])+y[ERKPP])
  v[23] = x[kf23]*y[RP]*y[PI3K] - x[kr23]*y[R_PI3K]
  v[24] = x[kf24]*y[R_PI3K] - x[kr24]*y[R_PI3K_act]
  v[25] = x[kf25]*y[R_PI3K_act] - x[kr25]*y[RP]*y[PI3K_act]
  v[26] = x[V26]*y[PI3K_act]/(x[K26]+y[PI3K_act])
  v[27] = x[kf27]*y[PI3K_act]*y[PI]/(x[K27]+y[PI])
  v[28] = x[V28]*y[PIP3]/(x[K28]+y[PIP3])
  v[29] = x[kf29]*y[PIP3]*y[Akt] - x[kr29]*y[Akt_PIP3]
  v[30] = x[V30]*y[Akt_PIP3]/(x[K30]*(1+y[Akt_PI_P]/x[K32])+y[Akt_PIP3])
  v[31] = x[kf31]*y[PP2A]*y[Akt_PI_P]/(x[K31]*(1+y[MEKP]/x[K16]+y[MEKPP]/x[K18]+y[Akt_PI_PP]/x[K33])+y[Akt_PI_P])
  v[32] = x[V32]*y[Akt_PI_P]/(x[K32]*(1+y[Akt_PIP3]/x[K30])+y[Akt_PI_P])
  v[33] = x[kf33]*y[PP2A]*y[Akt_PI_PP]/(x[K33]*(1+y[MEKP]/x[K16]+y[MEKPP]/x[K18]+y[Akt_PI_P]/x[K31])+y[Akt_PI_PP])
  v[34] = x[kf34]*y[RP] - x[kr34]*y[internalization]

  dydt = [0]*len(F_V)

  dydt[Akt] = -v[29]
  dydt[Akt_PIP3] = v[29] -v[30] + v[31]
  dydt[Akt_PI_P] = v[30] - v[31] - v[32] + v[33]
  dydt[Akt_PI_PP] = v[32] - v[33]
  dydt[ERK] = -v[19] + v[20]
  dydt[ERKP] = v[19] - v[20] - v[21] + v[22]
  dydt[ERKPP] = v[21] - v[22]
  dydt[GS] = -v[7] + v[9]
  dydt[HRG] = -v[1]
  dydt[internalization] = v[34]
  dydt[MEK] = -v[15] + v[16]
  dydt[MEKP] = v[15] - v[16] - v[17] + v[18]
  dydt[MEKPP] = v[17] - v[18]
  dydt[PI] = -v[27] + v[28]
  dydt[PI3K] = -v[23] + v[26]
  dydt[PI3K_act] = v[25] - v[26]
  dydt[PIP3] = v[27] - v[28] - v[29]
  dydt[R] = -v[1]
  dydt[RP] = v[3] - v[4] - v[5] + v[8] - v[23] + v[25] - v[34]
  dydt[R_HRG] = v[1] - 2*v[2]
  dydt[R_HRG2] = v[2] - v[3] + v[4]
  dydt[R_PI3K] = v[23] - v[24]
  dydt[R_PI3K_act] = v[24] - v[25]
  dydt[R_ShGS] = v[7] - v[8]
  dydt[R_ShP] = v[6] - v[7]
  dydt[R_Shc] = v[5] - v[6]
  dydt[Raf] = -v[13] + v[14]
  dydt[Raf_act] = v[13] - v[14]
  dydt[RasGDP] = -v[11] + v[12]
  dydt[RasGTP] = v[11] - v[12]
  dydt[ShGS] = v[8] - v[9]
  dydt[ShP] = v[9] - v[10]
  dydt[Shc] = -v[5] + v[10]

  return dydt
