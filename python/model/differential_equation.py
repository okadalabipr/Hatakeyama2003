from .name2idx import parameters as C
from .name2idx import variables as V

def diffeq(y,t,*x):

  v = [0]*35 #Rate equations

  v[1] = x[C.kf1]*y[V.R]*y[V.HRG] - x[C.kr1]*y[V.R_HRG]
  v[2] = x[C.kf2]*y[V.R_HRG]**2 - x[C.kr2]*y[V.R_HRG2]
  v[3] = x[C.kf3]*y[V.R_HRG2] - x[C.kr3]*y[V.RP]
  v[4] = x[C.V4]*y[V.RP]/(x[C.K4]+y[V.RP])
  v[5] = x[C.kf5]*y[V.RP]*y[V.Shc] - x[C.kr5]*y[V.R_Shc]
  v[6] = x[C.kf6]*y[V.R_Shc] - x[C.kr6]*y[V.R_ShP]
  v[7] = x[C.kf7]*y[V.R_ShP]*y[V.GS] - x[C.kr7]*y[V.R_ShGS]
  v[8] = x[C.kf8]*y[V.R_ShGS] - x[C.kr8]*y[V.ShGS]*y[V.RP]
  v[9] = x[C.kf9]*y[V.ShGS] - x[C.kr9]*y[V.GS]*y[V.ShP]
  v[10] = x[C.V10]*y[V.ShP]/(x[C.K10]+y[V.ShP])
  v[11] = x[C.kf11]*y[V.ShGS]*y[V.RasGDP]/(x[C.K11]+y[V.RasGDP])
  v[12] = x[C.V12]*y[V.RasGTP]/(x[C.K12]+y[V.RasGTP])
  v[13] = x[C.kf13]*y[V.RasGTP]*y[V.Raf]/(x[C.K13]+y[V.Raf])
  v[14] = x[C.kf14]*(y[V.Akt_PI_PP]+y[V.E])*y[V.Raf_act]/(x[C.K14]+y[V.Raf_act])
  v[15] = x[C.kf15]*y[V.Raf_act]*y[V.MEK]/(x[C.K15]*(1+y[V.MEKP]/x[C.K17])+y[V.MEK])
  v[16] = x[C.kf16]*y[V.PP2A]*y[V.MEKP]/(x[C.K16]*(1+y[V.MEKPP]/x[C.K18]+y[V.Akt_PI_P]/x[C.K31]+y[V.Akt_PI_PP]/x[C.K33])+y[V.MEKP])
  v[17] = x[C.kf17]*y[V.Raf_act]*y[V.MEKP]/(x[C.K17]*(1+y[V.MEK]/x[C.K15])+y[V.MEKP])
  v[18] = x[C.kf18]*y[V.PP2A]*y[V.MEKPP]/(x[C.K18]*(1+y[V.MEKP]/x[C.K16]+y[V.Akt_PI_P]/x[C.K31]+y[V.Akt_PI_PP]/x[C.K33])+y[V.MEKPP])
  v[19] = x[C.kf19]*y[V.MEKPP]*y[V.ERK]/(x[C.K19]*(1+y[V.ERKP]/x[C.K21])+y[V.ERK])
  v[20] = x[C.kf20]*y[V.MKP3]*y[V.ERKP]/(x[C.K20]*(1+y[V.ERKPP]/x[C.K22])+y[V.ERKP])
  v[21] = x[C.kf21]*y[V.MEKPP]*y[V.ERKP]/(x[C.K21]*(1+y[V.ERK]/x[C.K19])+y[V.ERKP])
  v[22] = x[C.kf22]*y[V.MKP3]*y[V.ERKPP]/(x[C.K22]*(1+y[V.ERKP]/x[C.K20])+y[V.ERKPP])
  v[23] = x[C.kf23]*y[V.RP]*y[V.PI3K] - x[C.kr23]*y[V.R_PI3K]
  v[24] = x[C.kf24]*y[V.R_PI3K] - x[C.kr24]*y[V.R_PI3K_act]
  v[25] = x[C.kf25]*y[V.R_PI3K_act] - x[C.kr25]*y[V.RP]*y[V.PI3K_act]
  v[26] = x[C.V26]*y[V.PI3K_act]/(x[C.K26]+y[V.PI3K_act])
  v[27] = x[C.kf27]*y[V.PI3K_act]*y[V.PI]/(x[C.K27]+y[V.PI])
  v[28] = x[C.V28]*y[V.PIP3]/(x[C.K28]+y[V.PIP3])
  v[29] = x[C.kf29]*y[V.PIP3]*y[V.Akt] - x[C.kr29]*y[V.Akt_PIP3]
  v[30] = x[C.V30]*y[V.Akt_PIP3]/(x[C.K30]*(1+y[V.Akt_PI_P]/x[C.K32])+y[V.Akt_PIP3])
  v[31] = x[C.kf31]*y[V.PP2A]*y[V.Akt_PI_P]/(x[C.K31]*(1+y[V.MEKP]/x[C.K16]+y[V.MEKPP]/x[C.K18]+y[V.Akt_PI_PP]/x[C.K33])+y[V.Akt_PI_P])
  v[32] = x[C.V32]*y[V.Akt_PI_P]/(x[C.K32]*(1+y[V.Akt_PIP3]/x[C.K30])+y[V.Akt_PI_P])
  v[33] = x[C.kf33]*y[V.PP2A]*y[V.Akt_PI_PP]/(x[C.K33]*(1+y[V.MEKP]/x[C.K16]+y[V.MEKPP]/x[C.K18]+y[V.Akt_PI_P]/x[C.K31])+y[V.Akt_PI_PP])
  v[34] = x[C.kf34]*y[V.RP] - x[C.kr34]*y[V.internalization]

  dydt = [0]*V.len_f_vars

  dydt[V.Akt] = -v[29]
  dydt[V.Akt_PIP3] = v[29] -v[30] + v[31]
  dydt[V.Akt_PI_P] = v[30] - v[31] - v[32] + v[33]
  dydt[V.Akt_PI_PP] = v[32] - v[33]
  dydt[V.ERK] = -v[19] + v[20]
  dydt[V.ERKP] = v[19] - v[20] - v[21] + v[22]
  dydt[V.ERKPP] = v[21] - v[22]
  dydt[V.GS] = -v[7] + v[9]
  dydt[V.HRG] = -v[1]
  dydt[V.internalization] = v[34]
  dydt[V.MEK] = -v[15] + v[16]
  dydt[V.MEKP] = v[15] - v[16] - v[17] + v[18]
  dydt[V.MEKPP] = v[17] - v[18]
  dydt[V.PI] = -v[27] + v[28]
  dydt[V.PI3K] = -v[23] + v[26]
  dydt[V.PI3K_act] = v[25] - v[26]
  dydt[V.PIP3] = v[27] - v[28] - v[29]
  dydt[V.R] = -v[1]
  dydt[V.RP] = v[3] - v[4] - v[5] + v[8] - v[23] + v[25] - v[34]
  dydt[V.R_HRG] = v[1] - 2*v[2]
  dydt[V.R_HRG2] = v[2] - v[3] + v[4]
  dydt[V.R_PI3K] = v[23] - v[24]
  dydt[V.R_PI3K_act] = v[24] - v[25]
  dydt[V.R_ShGS] = v[7] - v[8]
  dydt[V.R_ShP] = v[6] - v[7]
  dydt[V.R_Shc] = v[5] - v[6]
  dydt[V.Raf] = -v[13] + v[14]
  dydt[V.Raf_act] = v[13] - v[14]
  dydt[V.RasGDP] = -v[11] + v[12]
  dydt[V.RasGTP] = v[11] - v[12]
  dydt[V.ShGS] = v[8] - v[9]
  dydt[V.ShP] = v[9] - v[10]
  dydt[V.Shc] = -v[5] + v[10]

  return dydt
