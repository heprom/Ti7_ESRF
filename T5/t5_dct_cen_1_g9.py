import numpy as np

# data from the DCT experiment
# >> g9 = load('4_grains/phase_01/grain_0009.mat')
# >> edf_write(g9.proj.stack, 'g9_proj_stack.edf')
g_rod = [-0.4147, 0.1857, -0.0853]
g_center = [-0.2595, 0.0195, -0.1140]  # mm
g_hklplanes_exp = [
( 0,   1,  -1,   0),
( 0,  -1,   1,   0),
( 0,   1,  -1,   0),
( 0,  -1,   1,   0),
(-1,   1,   0,   0),
( 1,  -1,   0,   0),
( 1,  -1,   0,   0),
(-1,   1,   0,   0),
(-1,   0,   1,   0),
( 1,   0,  -1,   0),
( 1,   0,  -1,   0),
(-1,   0,   1,   0),
( 0,   0,   0,  -2),
( 0,   0,   0,   2),
( 0,   0,   0,   2),
( 0,  -1,   1,   1),
( 0,   1,  -1,  -1),
( 0,   1,  -1,  -1),
( 0,  -1,   1,   1),
( 0,  -1,   1,  -1),
(-1,   1,   0,   1),
( 1,  -1,   0,  -1),
( 1,  -1,   0,  -1),
(-1,   1,   0,   1),
(-1,   1,   0,  -1),
( 1,  -1,   0,   1),
( 1,  -1,   0,   1),
(-1,   1,   0,  -1),
(-1,   0,   1,   1),
( 1,   0,  -1,  -1),
( 1,   0,  -1,  -1),
(-1,   0,   1,   1),
(-1,   0,   1,  -1),
( 1,   0,  -1,   1),
( 1,   0,  -1,   1),
(-1,   0,   1,  -1),
( 0,  -1,   1,   2),
( 0,  -1,   1,   2),
( 1,  -1,   0,  -2),
( 1,  -1,   0,  -2),
( 1,  -1,   0,   2),
( 1,  -1,   0,   2),
(-1,   0,   1,   2),
( 1,   0,  -1,  -2),
(-1,   0,   1,   2),
( 1,   0,  -1,   2),
( 1,   0,  -1,   2),
( 1,  -2,   1,   0),
( 1,  -2,   1,   0),
(-1,  -1,   2,   0),
( 1,   1,  -2,   0),
(-1,  -1,   2,   0),
( 2,  -1,  -1,   0),
( 1,  -1,   0,   3),
( 1,  -1,   0,   3),
( 0,  -2,   2,   0),
( 0,  -2,   2,   0),
( 2,   0,  -2,   0),
( 2,   0,  -2,   0)]

# >> for i =1:54;disp(sprintf('(%d, %d, %d, %d),', g9.allblobs.hkl(g9.proj.ondet(g9.proj.included(i)), :)));end   
g_hklfamily = [
(0, 1, -1, 0),
(0, 1, -1, 0),
(0, 1, -1, 0),
(0, 1, -1, 0),
(0, 1, -1, 0),
(0, 1, -1, 0),
(0, 1, -1, 0),
(0, 1, -1, 0),
(0, 1, -1, 0),
(0, 1, -1, 0),
(0, 1, -1, 0),
(0, 1, -1, 0),
(0, 0, 0, 2),
(0, 0, 0, 2),
(0, 0, 0, 2),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 1),
(0, 1, -1, 2),
(0, 1, -1, 2),
(0, 1, -1, 2),
(0, 1, -1, 2),
(0, 1, -1, 2),
(0, 1, -1, 2),
(0, 1, -1, 2),
(0, 1, -1, 2),
(0, 1, -1, 2),
(0, 1, -1, 2),
(0, 1, -1, 2),
(-1, 2, -1, 0),
(-1, 2, -1, 0),
(-1, 2, -1, 0),
(-1, 2, -1, 0),
(-1, 2, -1, 0),
(-1, 2, -1, 0),
(0, 1, -1, 3),
(0, 1, -1, 3),
(0, 2, -2, 0),
(0, 2, -2, 0),
(0, 2, -2, 0),
(0, 2, -2, 0)]

g_om_exp = np.array([
    6.3377,
  186.3377,
  175.8579,
  355.8579,
  122.3998,
  302.3998,
  129.8959,
  309.8959,
   77.0679,
  257.0679,
   86.0415,
  266.0415,
   23.7080,
  203.7080,
   34.3280,
  164.2203,
  344.2203,
  139.9852,
  319.9852,
    6.8442,
  143.0998,
  323.0998,
  151.6375,
  331.6375,
   98.2002,
  278.2002,
  107.5237,
  287.5237,
  103.7008,
  283.7008,
  118.2228,
  298.2228,
   60.4991,
  240.4991,
   69.0070,
  249.0070,
  110.2224,
  247.1943,
  336.7919,
  168.5984,
  257.1751,
   90.5133,
  133.7021,
  155.8944,
  335.8944,
  229.9185,
   60.8646,
  321.2696,
  155.9578,
   35.5807,
  215.5807,
  234.5542,
  112.6246,
  241.6072,
   79.9200,
  191.6179,
  350.5857,
  252.5532,
   90.5561])

g_selected = np.array([
    0,
     1,
     1,
     1,
     1,
     1,
     1,
     1,
     0,
     1,
     0,
     0,
     1,
     1,
     1,
     1,
     1,
     1,
     1,
     0,
     1,
     1,
     1,
     1,
     1,
     1,
     1,
     1,
     1,
     1,
     1,
     0,
     1,
     1,
     1,
     1,
     1,
     0,
     0,
     0,
     1,
     0,
     1,
     0,
     1,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
     1,
     0,
     0,
     0,
     0,
     0])

#for om in g_om_exp:
#  print('full%04d.edf' % int(om/0.1))
# 44/44 spot in g4_proj_stack correspond to 109.2 deg -> full1092.edf (spot in u=1728 v=985 in image)

g_uv_exp = 1.0e+03 * np.array([
    [1.4601,    1.5926],
    [1.4621,    0.6621],
    [0.5479,    1.5516],
    [0.5765,    0.6301],
    [0.2220,    1.1935],
    [0.5272,    1.0124],
    [1.4963,    1.0155],
    [1.8279,    1.1966],
    [0.2619,    0.7211],
    [0.6658,    1.4727],
    [1.3686,    1.4955],
    [1.7982,    0.7369],
    [0.4026,    1.5885],
    [0.5756,    0.6572],
    [1.4626,    0.6217],
    [1.1975,    0.4353],
    [1.3431,    1.8262],
    [0.6393,    1.7851],
    [0.8954,    0.3949],
    [0.3038,    0.8624],
    [0.1861,    0.9671],
    [0.3972,    1.2497],
    [1.6234,    1.2399],
    [1.8696,    0.9552],
    [0.1653,    1.4267],
    [0.5351,    0.7783],
    [1.4962,    0.7817],
    [1.8741,    1.4303],
    [0.4170,    0.5306],
    [0.7687,    1.7127],
    [1.2694,    1.7030],
    [1.5463,    0.4826],
    [0.1186,    0.9662],
    [0.4620,    1.2436],
    [1.5878,    1.2466],
    [1.9045,    0.9688],
    [1.0703,    0.1834],
    [0.9535,    0.1893],
    [0.1910,    1.4950],
    [1.8270,    1.4699],
    [0.4167,    0.5501],
    [1.6349,    0.5302],
    [0.4180,    0.2848],
    [1.3809,    1.9217],
    [1.5943,    0.2361],
    [0.2306,    1.0075],
    [1.8451,    1.0087],
    [0.1356,    0.5313],
    [1.8908,    0.5592],
    [0.1349,    0.2438],
    [0.3760,    1.9319],
    [1.8865,    0.3206],
    [1.9538,    1.3982],
    [0.2247,    0.3185],
    [1.8177,    0.2899],
    [1.9218,    0.1999],
    [0.1196,    0.1429],
    [0.1234,    1.8631],
    [1.9301,    1.8671]])
