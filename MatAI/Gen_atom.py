# -*- coding: utf-8 -*-

symbol   =     ['X',  'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  
               'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  
               'Ca', 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 
               'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  
               'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 
               'Sn', 'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 
               'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 
               'Yb', 'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 
               'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 
               'Th', 'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 
               'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 
               'Ds', 'Rg', 'Cn', 'Uut','Fl', 'Uup','Lv', 'Uus']




period = [0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 
          2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
          4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 
          5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
          5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 
          6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
          6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
          6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 
          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
          7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
          7, 7, 7, 7, 7, 7, 7, 7, 7]




group = [0, 1, 18, 1, 2, 13, 14, 15, 16, 17, 
              18, 1, 2, 13, 14, 15, 16, 17, 18, 1, 
               2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
              12, 13, 14, 15, 16, 17, 18, 1, 2, 3, 
               4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
              14, 15, 16, 17, 18, 1, 2, 3, 3, 3, 
               3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
               3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
              12, 13, 14, 15, 16, 17, 18, 1, 2, 3, 
               3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
               3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 
              10, 11, 12, 13, 14, 15, 16, 17, 18]



mendeleevNumber = [0, 103, 1, 12, 77, 86, 95, 100, 101, 102, 
                   2, 11, 73, 80, 85, 90, 94, 99, 3, 10, 
                  16, 19, 51, 54, 57, 60, 61, 64, 67, 72, 
                  76, 81, 84, 89, 93, 98, 4, 9, 15, 25, 
                  49, 53, 56, 59, 62, 65, 69, 71, 75, 79, 
                  83, 88, 92, 97, 5, 8, 14, 33, 32, 31, 
                  30, 29, 28, 18, 27, 26, 24, 23, 22, 21, 
                  17, 20, 50, 52, 55, 58, 63, 66, 68, 70, 
                  74, 78, 82, 87, 91, 96, 6, 7, 13, 48, 
                  47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 
                  37, 36, 35, 34]




valence = [0, 1, 0, 1, 2, 3, 4, -3, -2, -1, 
           0, 1, 2, 3, 4, -3, -2, -1, 0, 1, 
           2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 
           2, 3, 4, 3, -2, -1, 0, 1, 2, 3, 
           4, 3, 4, 4, 3, 3, 2, 1, 2, 3, 
           4, 3, -2, -1, 0, 1, 2, 3, 3, 3, 
           3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
           3, 3, 4, 3, 4, 4, 4, 3, 2, 1, 
           2, 3, 4, 3, 2, -1, 0, 1, 2, 3, 
           3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
           3, 3, 3, 3]




atomicNum = {}
for an, symb in enumerate(symbol):
    atomicNum[symb] = an




radii = [0, 0.5292, 0.3113, 1.6282, 1.0855, 0.8141, 0.6513, 0.5427, 0.4652, 0.4071, 
               0.3618, 2.1649, 1.6711, 1.3607, 1.1476, 0.9922, 0.8738, 0.7807, 0.7056, 3.5598, 
               2.7479, 2.6106, 2.4861, 2.3732, 2.2701, 2.1754, 2.0885, 2.008, 1.9337, 1.8648, 
               1.8004, 1.5663, 1.3862, 1.2431, 1.1269, 1.0305, 0.9493, 4.8106, 3.7135, 3.5278, 
               3.3598, 3.2071, 3.0677, 2.9398, 2.8222, 2.7137, 2.6132, 2.5199, 2.433, 2.1167, 
               1.8732, 1.6799, 1.5228, 1.3926, 1.2828, 6.0615, 4.6788, 3.8102, 3.2133, 2.778, 
               2.4468, 2.1861, 1.9756, 1.802, 1.6565, 1.5328, 1.4262, 1.3335, 1.2521, 1.1801, 
               1.1159, 1.0583, 1.0079, 0.9594, 0.9165, 0.8773, 0.8413, 0.8182, 0.7776, 0.7492, 
               3.0636, 2.667, 2.3603, 2.1167, 1.9187, 1.7546, 1.6164, 7.2404, 5.5887, 5.3091, 
               5.0569, 3.7042, 3.2177, 2.8443, 2.3596, 2.1525, 2.1097, 1.8308, 1.7035, 1.5928, 
               1.4956, 1.4096, 1.3329, 1.3164]




radii2 = [0, 0.6216, 0.2908, 2.0872, 1.1329, 1.3216, 0.9422, 0.7032, 0.7673, 0.5734, 
                 0.4399, 2.2292, 1.4563, 1.9191, 1.3754, 1.0422, 1.0612, 0.8238, 0.6569, 2.6739, 
                 1.8637, 1.7336, 1.6586, 1.6829, 1.6789, 1.5151, 1.4252, 1.4258, 1.4774, 1.4601, 
                 1.1759, 1.9304, 1.4363, 1.1326, 1.1436, 0.9233, 0.7612, 2.7857, 2.0122, 1.7812, 
                 1.6532, 1.6454, 1.5922, 1.5509, 1.5322, 1.5138, 1.3398, 1.4917, 1.2349, 2.0063, 
                 1.555, 1.3046, 1.2487, 1.0606, 0.8979, 3.0019, 2.2149, 2.0707, 2.1243, 2.1494, 
                 2.1256, 2.1028, 2.0752, 2.064, 1.8968, 1.9985, 1.9728, 1.9424, 1.9154, 1.8891, 
                 1.8675, 2.172, 1.7507, 1.4578, 1.4406, 1.4609, 1.3111, 1.2474, 1.263, 1.2335, 
                 1.0764, 1.8934, 1.5388, 1.5712, 1.345, 1.1606, 1.0294, 2.805, 2.1846, 2.2371, 
                 1.8803, 1.9623, 1.8658, 1.8446, 1.919, 1.9454, 1.9356, 1.8693, 1.8486, 1.8131, 
                 1.7904, 1.7682, 1.7492, 2.5822]




covalentRadii = [0, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 
                  0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 
                  1.76, 1.7, 1.6, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 
                  1.22, 1.22, 1.2, 1.19, 1.2, 1.2, 1.16, 2.2, 1.95, 1.9, 
                  1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 
                  1.39, 1.39, 1.38, 1.39, 1.4, 2.44, 2.15, 2.07, 2.04, 2.03, 
                  2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.9, 
                  1.87, 1.87, 1.75, 1.7, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 
                  1.32, 1.45, 1.46, 1.48, 1.4, 1.5, 1.5, 2.6, 2.21, 2.15, 
                  2.06, 2, 1.96, 1.9, 1.87, 1.8, 1.69, 1.6, 1.6, 1.6, 
                  1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 
                  1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6]




vdwRadii = [0, 1.1, 1.4, 1.81, 1.53, 1.92, 1.7, 1.55, 1.52, 1.47, 
         1.54, 2.27, 1.73, 1.84, 2.1, 1.8, 1.8, 1.75, 1.88, 2.75, 
         2.31, 2.3, 2.15, 2.05, 2.05, 2.05, 2.05, 2, 2, 2, 
          2.1, 1.87, 2.11, 1.85, 1.9, 1.83, 2.02, 3.03, 2.49, 2.4, 
          2.3, 2.15, 2.1, 2.05, 2.05, 2, 2.05, 2.1, 2.2, 2.2, 
         1.93, 2.17, 2.06, 1.98, 2.16, 3.43, 2.68, 2.5, 2.48, 2.47, 
         2.45, 2.43, 2.42, 2.4, 2.38, 2.37, 2.35, 2.33, 2.32, 2.3, 
         2.28, 2.27, 2.25, 2.2, 2.1, 2.05, 2, 2, 2.05, 2.1, 
         2.05, 1.96, 2.02, 2.07, 1.97, 2.02, 2.2, 3.48, 2.83, 2, 
          2.4, 2, 2.3, 2, 2, 2, 2, 2, 2, 2, 
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
            2, 2, 2, 2, 2, 2, 2, 2]



mass = [0, 1.00794, 4.0026, 6.941, 9.01218, 10.811, 12.0107, 14.0067, 15.9994, 18.9984, 
            20.1797, 22.9898, 24.305, 26.9815, 28.0855, 30.9738, 32.065, 35.453, 39.948, 39.0983, 
            40.078, 44.9559, 47.867, 50.9415, 51.9961, 54.938, 55.845, 58.9332, 58.6934, 63.546, 
            65.409, 69.723, 72.64, 74.9216, 78.96, 79.904, 83.798, 85.4678, 87.62, 88.9059, 
            91.224, 92.9064, 95.94, 98, 101.07, 102.906, 106.42, 107.868, 112.411, 114.818, 
            118.71, 121.76, 127.6, 126.904, 131.293, 132.905, 137.327, 138.905, 140.116, 140.908, 
            144.242, 145, 150.36, 151.964, 157.25, 158.925, 162.5, 164.93, 167.259, 168.934, 
            173.04, 174.967, 178.49, 180.948, 183.84, 186.207, 190.23, 192.217, 195.084, 196.967, 
            200.59, 204.383, 207.2, 208.98, 210, 210, 220, 223, 226, 227, 
            232.038, 231.036, 238.029, 237, 244, 243, 247, 247, 251, 252, 
            257, 258, 259, 262, 265, 268, 271, 270, 277, 276, 
            281, 280, 285, 284, 289, 288, 293, 294]



electronAffinity = [0, 0.754195, -0.52, 0.618049, -0.52, 0.279723, 1.26212, -0.000725, 1.46111, 3.40119, 
                -1.2, 0.547926, -0.415, 0.43283, 1.38952, 0.746607, 2.0771, 3.61272, -1, 0.501459, 
             0.02455, 0.188, 0.084, 0.52766, 0.67584, -0.52, 0.153236, 0.66226, 1.15716, 1.23578, 
               -0.62, 0.43, 1.23268, 0.8048, 2.0206, 3.36359, -0.62, 0.485916, 0.05206, 0.307, 
              0.4333, 0.9174, 0.7473, 0.55, 1.04638, 1.14289, 0.56214, 1.30447, -0.725, 0.3, 
             1.11207, 1.0474, 1.97087, 3.05905, -0.83, 0.47163, 0.14462, 0.47, 0.65, 0.962, 
               1.916, 0.129, 0.162, 0.864, 0.137, 1.165, 0.352, 0.338, 0.312, 1.029, 
               -0.02, 0.346, 0.017, 0.323, 0.81626, 0.060396, 1.1, 1.56436, 2.1251, 2.30861, 
               -0.52, 0.377, 0.356743, 0.942362, 1.9, 2.3, -0.725, 0.486, 0.1, 0.35, 
                1.17, 0.55, 0.53, 0.48, -0.5, 0.1, 0.28, -1.72, -1.01, -0.3, 
                0.35, 0.98, -2.33, -0.31]



electronNegativity = [0, 2.2, 4.16, 0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 
                 4.787, 0.93, 1.31, 1.61, 1.9, 2.19, 2.58, 3.16, 3.242, 0.82, 
                     1, 1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.9, 
                  1.65, 1.81, 2.01, 2.18, 2.55, 2.96, 3, 0.82, 0.95, 1.22, 
                  1.33, 1.6, 2.16, 1.9, 2.2, 2.28, 2.2, 1.93, 1.69, 1.78, 
                  1.96, 2.05, 2.1, 2.66, 2.6, 0.79, 0.89, 1.1, 1.12, 1.13, 
                  1.14, 1, 1.17, 1, 1.2, 1, 1.22, 1.23, 1.24, 1.25, 
                     1, 1.27, 1.3, 1.5, 2.36, 1.9, 2.2, 2.2, 2.28, 2.54, 
                     2, 1.62, 2.33, 2.02, 2, 2.2, 2.6, 0.67, 0.9, 1.1, 
                   1.3, 1.5, 1.38, 1.36, 1.28, 1.3, 1.3, 1.3, 1.3, 1.3, 
                   1.3, 1.3, 1.3]



electronNegativity2 = [0, 0.2638, 0.442712, 0.105039, 0.144986, 0.184886, 0.224776, 0.26493, 0.304575, 0.344443, 
                0.38439, 0.093214, 0.121644, 0.150078, 0.178503, 0.206931, 0.23596, 0.263803, 0.2922, 0.0981896, 
               0.115412, 0.119383, 0.123364, 0.127334, 0.131305, 0.135284, 0.139253, 0.143236, 0.147207, 0.151172, 
               0.155152, 0.172377, 0.189589, 0.206821, 0.224033, 0.241258, 0.258481, 0.104686, 0.118508, 0.121699, 
               0.124889, 0.173447, 0.177911, 0.182377, 0.186842, 0.19131, 0.195772, 0.08731, 0.09789, 0.10033, 
                0.10277, 0.12149, 0.13207, 0.128078, 0.131267, 0.134459, 0.137649, 0.140838, 0.144028, 0.147217, 
               0.150407, 0.16423, 0.178052, 0.191877, 0.205699, 0.219518, 0.233351, 0.154213, 0.158679, 0.163142, 
                0.16761, 0.17208, 0.176539, 0.181003, 0.185468, 0.189935, 0.1944, 0.198863, 0.20333, 0.207795, 
               0.212261, 0.216724, 0.22119, 0.22565, 0.229987, 0.234581, 0.23905, 0.243516, 0.247984, 0.25106, 
                0.25691, 0.26137, 0.169, 0.14265, 0.16137, 0.17194, 0.17439, 0.19315, 0.20369, 0.21427, 
                0.22485, 0.23542, 0.246, 0.24844]




ionizationEnergy = [0, 13.5984, 24.5874, 5.39172, 9.3227, 8.29803, 11.2603, 14.5341, 13.6181, 17.4228, 
                21.5646, 5.13908, 7.64624, 5.98577, 8.15169, 10.4867, 10.36, 12.9676, 15.7596, 4.34066, 
                6.11316, 6.5615, 6.8281, 6.7462, 6.7665, 7.43402, 7.9024, 7.881, 7.6398, 7.72638, 
                 9.3942, 5.9993, 7.8994, 9.7886, 9.75238, 11.8138, 13.9996, 4.17713, 5.6949, 6.2171, 
                 6.6339, 6.75885, 7.09243, 7.28, 7.3605, 7.4589, 8.3369, 7.5762, 8.9938, 5.78636, 
                 7.3439, 8.6084, 9.0096, 10.4513, 12.1298, 3.8939, 5.2117, 5.5769, 5.5387, 5.473, 
                  5.525, 5.582, 5.6436, 5.6704, 6.1501, 5.8638, 5.9389, 6.0215, 6.1077, 6.18431, 
                6.25416, 5.4259, 6.82507, 7.5496, 7.864, 7.8335, 8.4382, 8.967, 8.9587, 9.2255, 
                10.4375, 6.1082, 7.41666, 7.2856, 8.417, 9.31751, 10.7485, 4.0727, 5.2784, 5.17, 
                 6.3067, 5.89, 6.19405, 6.2657, 6.0262, 5.9738, 5.9915, 6.1979, 6.2817, 6.42, 
                    6.5, 6.58, 6.65, 4.9, 6.01]




bolingPoint = [0, 20.28, 4.22, 1615, 2742, 4200, 4300, 77.36, 90.2, 85.03, 
                 27.07, 1156, 1363, 2792, 3173, 550, 717.87, 239.11, 87.3, 1032, 
                  1757, 3103, 3560, 3680, 2944, 2334, 3134, 3200, 3186, 3200, 
                  1180, 2477, 3093, 887, 958, 332, 119.93, 961, 1655, 3609, 
                  4682, 5017, 4912, 4538, 4423, 3968, 3236, 2435, 1040, 2345, 
                  2875, 1860, 1261, 457.4, 165.1, 944, 2143, 3743, 3633, 3563, 
                  3373, 3273, 2076, 1800, 3523, 3503, 2840, 2993, 3141, 2223, 
                  1469, 3675, 4876, 5731, 5828, 5869, 5285, 4701, 4098, 3129, 
                629.88, 1746, 2022, 1837, 1235, 503, 211.3, 871, 2010, 3573, 
                  5093, 4300, 4200, 4273, 3503, 2880, 3383, 3000, 3000, 3000, 
                  3000, 3000, 3000, 3000]



meltingPoint = [0, 14.01, 0.95, 453.69, 1560, 2349, 3800, 63.05, 54.8, 53.53, 
                  24.56, 370.87, 923, 933.47, 1687, 317.3, 388.36, 171.6, 83.8, 336.53, 
                   1115, 1814, 1941, 2183, 2180, 1519, 1811, 1768, 1728, 1357.77, 
                 692.68, 302.91, 1211.4, 1090, 494, 265.8, 115.79, 312.46, 1050, 1799, 
                   2128, 2750, 2896, 2430, 2607, 2237, 1828.05, 1234.93, 594.22, 429.75, 
                 505.08, 903.78, 722.66, 386.85, 161.4, 301.59, 1000, 1193, 1068, 1208, 
                   1297, 1373, 1345, 1099, 1585, 1629, 1680, 1734, 1802, 1818, 
                   1097, 1925, 2506, 3290, 3695, 3459, 3306, 2739, 2041.4, 1337.33, 
                 234.32, 577, 600.61, 544.4, 527, 575, 202, 300, 973, 1323, 
                   2115, 1841, 1405.3, 910, 912.5, 1449, 1613, 1259, 1173, 1133, 
                   1800, 1100, 1100, 1900]




moleVolume = [0, 11.42, 21, 13.02, 4.85, 4.39, 5.29, 13.54, 17.36, 11.2, 
                13.23, 23.78, 14, 10, 12.06, 17.02, 15.53, 17.39, 22.56, 45.94, 
                 26.2, 15, 10.64, 8.32, 7.23, 7.35, 7.09, 6.67, 6.59, 7.11, 
                 9.16, 11.8, 13.63, 12.95, 16.42, 19.78, 27.99, 55.76, 33.94, 19.88, 
                14.02, 10.83, 9.38, 8.63, 8.17, 8.28, 8.56, 10.27, 13, 15.76, 
                16.29, 18.19, 20.46, 25.72, 35.92, 70.94, 38.16, 22.39, 20.69, 20.8, 
                20.59, 20.23, 19.98, 28.97, 19.9, 19.3, 19.01, 18.74, 18.46, 19.1, 
                24.84, 17.78, 13.44, 10.85, 9.47, 8.86, 8.42, 8.52, 9.09, 10.21, 
                14.09, 17.22, 18.26, 21.31, 22.97, 33, 50.5, 77, 41.09, 22.55, 
                 19.8, 15.18, 12.49, 11.59, 12.29, 17.63, 18.05, 16.84, 16.5, 28.52, 
                   12, 12, 12, 12, 16, 12, 12, 10, 10, 10, 
                   10, 12, 17, 18, 21, 22, 26, 41, 52]




thermConductivity = [0, 0.1805, 0.1513, 85, 190, 27, 140, 0.02583, 0.02658, 0.0277, 
                 0.0491, 140, 160, 235, 150, 0.236, 0.205, 0.0089, 0.01772, 100, 
                    200, 16, 22, 31, 94, 7.8, 80, 100, 91, 400, 
                    120, 29, 60, 50, 0.52, 0.12, 0.00943, 58, 35, 17, 
                     23, 54, 139, 51, 120, 150, 72, 430, 97, 82, 
                     67, 24, 3, 0.449, 0.00565, 36, 18, 13, 11, 13, 
                     17, 15, 13, 14, 11, 11, 11, 16, 15, 17, 
                     39, 16, 23, 57, 170, 48, 88, 150, 72, 320, 
                    8.3, 46, 35, 8, 20, 2, 0.00361, 15, 19, 12, 
                     54, 47, 27, 6, 6, 10, 8.8, 10, 10, 10, 
                     10, 10, 10, 10, 23, 58, 19]




orbitalExponent = [0, 1, 1.7, 0.65, 0.975, 1.3, 1.625, 1.95, 2.275, 2.6, 
                        2.925, 0.7333, 0.95, 1.1667, 1.3833, 1.6, 1.8167, 2.0333, 2.25, 0.5946, 
                       0.7703, 0.8108, 0.8514, 0.8919, 0.9324, 0.973, 1.0135, 1.0541, 1.0946, 1.1351, 
                       1.1757, 1.3514, 1.527, 1.7027, 1.8784, 2.0541, 2.2297, 0.55, 0.7125, 0.75, 
                       0.7875, 0.825, 0.8625, 0.9, 0.9375, 0.975, 1.0125, 1.05, 1.0875, 1.25, 
                       1.4125, 1.575, 1.7375, 1.9, 2.0625, 0.5238, 0.6786, 0.8333, 0.9881, 1.1429, 
                       1.2976, 1.4524, 1.6071, 1.7619, 1.9167, 2.0714, 2.2262, 2.381, 2.5357, 2.6905, 
                       2.8452, 3, 3.15, 3.3095, 3.4643, 3.619, 3.7738, 3.9286, 4.0833, 4.2381, 
                       4.3929, 1.1905, 1.3452, 1.5, 1.6548, 1.8095, 1.9643, 0.5116, 0.6628, 0.6977, 
                       0.7326, 1, 1.1512, 1.3023, 1.5698, 1.7209, 1.7558, 2.0233, 2.1744, 2.3256, 
                       2.4767, 2.6279, 2.7791, 2.814]




polarizability = [0, 0.6669, 0.1358, 19.4239, 5.7558, 2.428, 1.2432, 0.7193, 0.453, 0.3034, 
                   0.2131, 45.659, 21, 11.337, 6.8011, 4.3955, 3.0023, 2.1412, 1.5808, 202.997, 
                  93.3717, 80.0633, 69.1712, 60.1472, 52.6438, 46.3265, 40.9939, 36.4555, 32.5372, 29.1769, 
                  26.2615, 17.2917, 11.9864, 8.6443, 6.4397, 4.9244, 3.8497, 500.968, 230.446, 197.571, 
                  170.668, 148.44, 129.913, 114.332, 101.152, 89.9286, 80.3028, 72.005, 64.8095, 42.6767, 
                  29.5777, 21.3335, 15.8906, 12.1532, 9.4993, 1002.2, 460.91, 248.918, 149.288, 96.4738, 
                  65.9186, 47.0135, 34.6984, 26.3316, 20.4544, 16.2057, 13.0543, 10.6707, 8.8334, 7.3955, 
                    6.253, 5.3338, 4.6075, 3.9739, 3.4643, 3.0385, 2.6796, 2.3756, 2.1175, 1.8924, 
                  129.646, 85.3653, 59.1717, 42.6767, 31.7858, 24.3079, 19.0046, 1705, 785.498, 673.403, 
                  581.682, 228.716, 149.916, 103.547, 59.1266, 44.8789, 42.2547, 27.5918, 22.2453, 18.1843, 
                  15.0542, 12.6038, 10.6563, 10.2631]




globalHardness = [0, 13.588, 22.383, 4.4164, 6.6244, 8.8328, 11.0407, 13.25, 15.4574, 17.6634, 
             19.875, 3.3215, 4.303, 5.2846, 6.2659, 7.2473, 8.2293, 9.2107, 10.191, 2.02, 
             2.6162, 2.7545, 2.8924, 3.03, 3.1676, 3.3055, 3.443, 3.5811, 3.7187, 3.8561, 
              3.994, 4.5909, 5.1874, 5.7846, 6.381, 6.978, 7.5748, 1.4948, 1.9364, 2.0383, 
             2.1402, 2.2422, 2.344, 2.446, 2.5479, 2.6498, 2.7517, 2.8536, 2.9555, 3.3972, 
             3.8388, 4.2805, 4.7221, 5.0632, 5.6056, 1.1863, 1.5369, 1.8873, 2.2378, 2.5885, 
             2.9389, 3.2893, 3.6398, 3.9905, 4.341, 4.6913, 5.0419, 5.3924, 5.743, 6.0934, 
             6.4439, 6.7947, 7.1344, 7.4951, 7.8459, 8.1965, 8.5472, 8.8973, 9.2474, 9.598, 
             2.3456, 2.6962, 3.0466, 3.3972, 3.7477, 4.0983, 4.4487, 0.9913, 1.2867, 1.3544, 
             1.4222, 1.9413, 2.2348, 2.5281, 3.0473, 3.3407, 3.4084, 3.9277, 4.2212, 4.5146, 
              4.808, 5.1013, 5.3949, 5.4629]



electrophilicity = [0, 5.6326, 8.85722, 2.52302, 3.2717, 4.02042, 4.76896, 5.51718, 6.26638, 7.0143, 
              7.76411, 2.15186, 2.48467, 2.81733, 3.14998, 3.48289, 3.8155, 4.14814, 4.4809, 1.76613, 
              1.98489, 2.03533, 2.08586, 2.13632, 2.18671, 2.23725, 2.28773, 2.33827, 2.38868, 2.43909, 
              2.48967, 2.7084, 2.92701, 3.15499, 3.36456, 3.58342, 3.80184, 1.65924, 1.8464, 1.88958, 
              1.93276, 1.97596, 2.01913, 2.06233, 2.10553, 2.1487, 2.19188, 2.23509, 2.27826, 2.46546, 
               2.6526, 2.83973, 3.02688, 3.21404, 3.40125, 1.60033, 1.77014, 1.9398, 2.10961, 2.27942, 
              2.44907, 2.6189, 2.78858, 2.95836, 3.12818, 3.29786, 3.46763, 3.63739, 3.80723, 3.97692, 
              4.14654, 4.31628, 4.4809, 4.65582, 4.82555, 4.99571, 5.16488, 5.33536, 5.5048, 5.67477, 
              5.84382, 2.3316, 2.50127, 2.67104, 2.84094, 3.01059, 3.18039, 1.57391, 1.7359, 1.77328, 
              1.81068, 2.09713, 2.25913, 2.42105, 2.70759, 2.91347, 2.90691, 3.19342, 3.35539, 3.51729, 
              3.67917, 3.84127, 4.00323, 4.04079]




atomizationEnthalpy = [0, 218, 0, 159, 324, 563, 717, 473, 249, 79, 
                  0, 107, 146, 326, 456, 315, 279, 122, 0, 89, 
                178, 378, 471, 515, 397, 281, 415, 426, 431, 338, 
                131, 277, 377, 302, 227, 112, 0, 81, 164, 425, 
                605, 733, 659, 661, 652, 556, 377, 285, 112, 243, 
                302, 262, 197, 107, 0, 76, 182, 431, 423, 356, 
                328, 350, 207, 175, 398, 389, 290, 301, 317, 232, 
                152, 428, 621, 782, 860, 776, 789, 671, 565, 368, 
                 64, 182, 195, 207, 142, 0, 0, 64, 159, 406, 
                598, 607, 536]




fusionEnthalpy = [0, 0.558, 0.02, 3, 7.95, 50, 117, 0.36, 0.222, 0.26, 
               0.34, 2.6, 8.7, 10.7, 50.2, 0.64, 1.73, 3.2, 1.18, 2.33, 
               8.54, 16, 18.7, 22.8, 20.5, 13.2, 13.8, 16.2, 17.2, 13.1, 
               7.35, 5.59, 31.8, 27.7, 5.4, 5.8, 1.64, 2.19, 8, 11.4, 
                 21, 26.8, 36, 23, 25.7, 21.7, 16.7, 11.3, 6.3, 3.26, 
                  7, 19.7, 17.5, 7.76, 2.3, 2.09, 8, 6.2, 5.5, 6.9, 
                7.1, 7.7, 8.6, 9.2, 10, 10.8, 11.1, 17, 19.9, 16.8, 
                7.7, 22, 25.5, 36, 35, 33, 31, 26, 20, 12.5, 
               2.29, 4.2, 4.77, 10.9, 13, 6, 3, 2, 8, 14, 
                 16, 15, 14, 10, 2.8]



vapEnthalpy = [0, 0.452, 0.083, 147, 297, 507, 715, 2.79, 3.41, 3.27, 
            1.75, 97.7, 128, 293, 359, 12.4, 9.8, 10.2, 6.5, 76.9, 
             155, 318, 425, 453, 339, 220, 347, 375, 378, 300, 
             119, 256, 334, 32.4, 26, 14.8, 9.02, 72, 137, 380, 
             580, 690, 600, 550, 580, 495, 380, 255, 100, 230, 
             290, 68, 48, 20.9, 12.64, 65, 140, 400, 350, 330, 
             285, 290, 175, 175, 305, 295, 280, 265, 285, 250, 
             160, 415, 630, 735, 800, 705, 630, 560, 490, 330, 
            59.2, 165, 178, 160, 100, 40, 17, 65, 125, 400, 
             530, 470, 420, 335, 325]



bindingEnergy = [0, 13.6, 24.6, 54.7, 111.5, 188, 284.2, 409.9, 543.1, 696.7, 
             870.2, 1070.8, 1303, 1559, 1839, 2145.5, 2472, 2822, 3205.9, 3608.4, 
            4038.5, 4492, 4966, 5465, 5989, 6539, 7112, 7709, 8333, 8979, 
              9659, 10367, 11103, 11867, 12658, 13474, 14326, 15200, 16105, 17038, 
             17998, 18986, 20000, 21044, 22117, 23220, 24350, 25514, 26711, 27940, 
             29200, 30491, 31814, 33169, 34561, 35985, 37441, 38925, 40443, 41991, 
             43569, 45184, 46834, 48519, 50239, 51996, 53789, 55618, 57486, 59390, 
             61332, 63314, 65351, 67416, 69525, 71676, 73871, 76111, 78395, 80725, 
             83102, 85530, 88005, 90526, 93105, 95730, 98404, 101137, 103922, 106755, 
            109651, 112601, 115606]



density = [0, 88, 214, 535, 1848, 2460, 2267, 1026, 1495, 1700, 
        1444, 968, 1738, 2700, 2330, 1823, 1960, 2030, 1616, 856, 
        1550, 2985, 4507, 6110, 7140, 7470, 7874, 8900, 8908, 8920, 
        7140, 5904, 5323, 5727, 4819, 4050, 2155, 1532, 2630, 4472, 
        6511, 8570, 10280, 11500, 12370, 12450, 12023, 10490, 8650, 7310, 
        7310, 6697, 6240, 4940, 3640, 1879, 3510, 6146, 6689, 6640, 
        6800, 7264, 7353, 5244, 7901, 8219, 8551, 8795, 9066, 9321, 
        6570, 9841, 13310, 16650, 19250, 21020, 22610, 22650, 21090, 19300, 
       14190, 11850, 11340, 9780, 9196, 6400, 4400, 2900, 5000, 10070, 
       11724, 15370, 19050, 20450, 19816, 13780, 13510, 14780, 15100, 13500, 
        1400, 1400, 1400, 1400, 17000, 21600, 23200, 27200, 28600, 28200, 
       27400, 24400, 16800, 16000, 14000, 13000, 11200, 7200, 5700]




class atom:
    def __init__(self,symbol,num=1,x=0.0,y=0.0,z=0.0):
        self.name = symbol
        _an = atomicNum[symbol]
        self.atomicNum = _an
        self.num = num
        self.x = x
        self.y = y
        self.z = z
    def get_mendeleevNumber(self):
        return mendeleevNumber[self.atomicNum]
    def get_period(self):
        return period[self.atomicNum]
    def get_group(self):
        return group[self.atomicNum]
    def get_valence(self):
        return valence[self.atomicNum]
    def get_radii(self):
        return radii[self.atomicNum]
    def get_radii2(self):
        return radii2[self.atomicNum]
    def get_covalentRadii(self):
        return covalentRadii[self.atomicNum]
    def get_vdwRadii(self):
        return vdwRadii[self.atomicNum]
    def get_mass(self):
        return mass[self.atomicNum]
    def get_density(self):
        return density[self.atomicNum]
    def get_electronAffinity(self):
        return electronAffinity[self.atomicNum]
    def get_electronNegativity(self):
        return electronNegativity[self.atomicNum]
    def get_electronNegativity2(self):
        return electronNegativity2[self.atomicNum]
    def get_ionizationEnergy(self):
        return ionizationEnergy[self.atomicNum]
    def get_bolingPoint(self):
        return bolingPoint[self.atomicNum]
    def get_meltingPoint(self):
        return meltingPoint[self.atomicNum]
    def get_moleVolume(self):
        return moleVolume[self.atomicNum]
    def get_thermConductivity(self):
        return thermConductivity[self.atomicNum]
    def get_orbitalExponent(self):
        return orbitalExponent[self.atomicNum]
    def get_polarizability(self):
        return polarizability[self.atomicNum]
    def get_globalHardness(self):
        return globalHardness[self.atomicNum]
    def get_electrophilicity(self):
        return electrophilicity[self.atomicNum]
    def get_atomizationEnthalpy(self):
        return atomizationEnthalpy[self.atomicNum]
    def get_fusionEnthalpy(self):
        return fusionEnthalpy[self.atomicNum]
    def get_vapEnthalpy(self):
        return vapEnthalpy[self.atomicNum]
    def get_bindingEnergy(self):
        return bindingEnergy[self.atomicNum]       

    def get_property(self):
        ppty = [self.atomicNum,
                self.get_mendeleevNumber(),
                self.get_period(),
                self.get_group(),
                self.get_mass(),
                self.get_density(),
                self.get_valence(),
                self.get_radii2(),
                self.get_covalentRadii(),
                self.get_vdwRadii(),
                self.get_electronAffinity(),
                self.get_electronNegativity2(),
                self.get_ionizationEnergy(),
                self.get_bolingPoint(),
                self.get_meltingPoint(),
                self.get_moleVolume(),
                self.get_thermConductivity(),
                self.get_orbitalExponent(),
                self.get_polarizability(),
                self.get_globalHardness(),
                self.get_electrophilicity(),
                self.get_atomizationEnthalpy(),
                self.get_fusionEnthalpy(),
                self.get_vapEnthalpy(),
                self.get_bindingEnergy()]
        return ppty
                
                
if __name__ == '__main__':
    at = atom("Fe")
    print(at.get_radii())
    ats = [atom("Al"),atom("Cr")]
    print(ats[1].get_radii())


