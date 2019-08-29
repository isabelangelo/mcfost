
### DUST MASS WEIGHTS ###

#non-normalized (these have negatives)
# mass_lupus = [1.051021015953335,
#  0.766021015953335,
#  0.48102101595333513,
#  0.196021015953335,
#  -0.08897898404666515]
# mass_taurus = [1.059846683459652,
#  0.7598466834596516,
#  0.4598466834596511,
#  0.1598466834596508,
#  -0.14015331654034957]

#normalized
mass_lupus_norm = [0.42140561, 0.3071352 , 0.1928648 , 0.07859439, 0.        ]
mass_taurus_norm = [0.43447259, 0.31149086, 0.18850914, 0.06552741, 0.        ]

# mass with Mdust=1e-7 truncated
mass_taurus_truncated = mass_taurus_norm.copy()
mass_taurus_truncated[0] = 0

### H0 WIEIGHTS ###

#non-normalized
H0_hist = [6304, 6586, 1267,  324]
H0_gaussian = [1202.9117685443744, 1415.4914349640962, 273.25111495679084, 8.653620875746345]

#normalized
H0_hist_norm = [0.43532905, 0.45480285, 0.08749396, 0.02237415]
H0_gaussian_norm = [0.41475312, 0.48804867, 0.09421452, 0.00298369]


### Rc WEIGHTS ###
Rc_norm = [1/3., 1/3., 1/3., 0]

### f_exp WEIGHTS ###
f_exp_norm = [0, 1/3., 1/3., 1/3.]