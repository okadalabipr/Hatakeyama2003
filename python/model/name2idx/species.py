NAMES = [
    'Akt',
    'Akt_PIP3',
    'Akt_PI_P',
    'Akt_PI_PP',
    'ERK',
    'ERKP',
    'ERKPP',
    'GS',
    'HRG',
    'internalization',
    'MEK',
    'MEKP',
    'MEKPP',
    'PI',
    'PI3K',
    'PI3K_act',
    'PIP3',
    'R',
    'RP',
    'R_HRG',
    'R_HRG2',
    'R_PI3K',
    'R_PI3K_act',
    'R_ShGS',
    'R_ShP',
    'R_Shc',
    'Raf',
    'Raf_act',
    'RasGDP',
    'RasGTP',
    'ShGS',
    'ShP',
    'Shc',
    'E',
    'MKP3',
    'PP2A',
]

for idx, name in enumerate(NAMES):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

NUM = len(NAMES)