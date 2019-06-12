#ifndef D_Z_CONFIG_H
#define D_Z_CONFIG_H

// I am aware this isn't very human-readable -- python has issues with line sizes, so the lines had to be truncated in some fashion
const char * cosmos[5] =  {"PLANCK15", "PLANCK13", "WMAP9", "WMAP7", "WMAP5"};
const int num_cosmologies = 5;const int num_segments[5] = { 3, 3, 3, 3, 3};
const int interp_degree[5] = { 12, 12, 12, 12, 12};
const double boundaries_Z[5][4] = {{0.00027144176165949066
, 0.00027144176165949066
, 0.07368062997280773
, 20.000000000000004
},{0.00027144176165949066
, 0.00027144176165949066
, 0.07368062997280773
, 20.000000000000004
},{0.00027144176165949066
, 0.00027144176165949066
, 0.07368062997280773
, 20.000000000000004
},{0.00027144176165949066
, 0.00027144176165949066
, 0.07368062997280773
, 20.000000000000004
},{0.00027144176165949066
, 0.00027144176165949066
, 0.07368062997280773
, 20.000000000000004
}};
const double boundaries_D[5][4] = {{1.2015524178010442
, 1.2015524178010442
, 344.0498142197095
, 230421.91551573356
},{1.201020615326616
, 1.201020615326616
, 343.9047300972533
, 230427.0888269416
},{1.17417097464182
, 1.17417097464182
, 336.6234405707774
, 231610.64191786345
},{1.1561615382304566
, 1.1561615382304566
, 331.7253154763896
, 232528.33987307717
},{1.1594542653630548
, 1.1594542653630548
, 332.57833550529614
, 231611.2292384872
}};
const double COEFF_VEC_DZ[5][3][12] =  {{{ -2.4916278904528094e-11
, 1.0074682884645874e-09
, 0.00022593977849538305
, 1.464141132428875e-07
, -8.232937802852727e-07
, 2.6959557250839955e-06
, -6.136907535217196e-06
, 9.335538393977773e-06
, -9.376606674281895e-06
, 5.963656562234233e-06
, -2.1745207423992716e-06
, 3.461147455366076e-07
},{ 7.744302358729577e-09
, -1.8420646340129994e-08
, 0.00022597406280918379
, -9.208274973545955e-09
, -3.633242812566835e-08
, -5.830675867067343e-10
, 9.474078367339846e-11
, -6.960292055329224e-12
, 4.0203433116113996e-13
, -1.5360475121306272e-14
, 3.352769238341968e-16
, -3.1715447860533055e-18
},{ 0.004560877458856844
, -0.0010488155718547695
, 0.0003166704087496086
, -3.826659151626778e-06
, 3.78150930364314e-08
, -2.4890905790782446e-10
, 1.1173317121815323e-12
, -3.4338475416696076e-15
, 7.105409999825397e-18
, -9.457921645468802e-21
, 7.309765914691665e-24
, -2.4916848759769718e-27
}},{{ 3.0608026655942855e-11
, -1.2254593538442233e-09
, 0.0002260762839393041
, -1.741490174832977e-07
, 8.833437684210174e-07
, -3.141494921386785e-06
, 7.089610382673379e-06
, -1.0703752232593744e-05
, 1.0680840573656513e-05
, -6.754913267653557e-06
, 2.4510380629342868e-06
, -3.884812704553851e-07
},{ -6.844420213734554e-09
, 1.637927871382273e-08
, 0.00022604052012505408
, 8.277764015710223e-09
, -4.1884988414153865e-08
, 5.296346038792078e-10
, -5.428903092858041e-11
, 6.3943178919544994e-12
, -3.885245871052655e-13
, 1.4290171936510378e-14
, -3.0259366344120127e-16
, 2.821055941325703e-18
},{ 0.0045453019264830846
, -0.001047425896196084
, 0.00031675125387406315
, -3.828354388678864e-06
, 3.783108685749574e-08
, -2.489990937081228e-10
, 1.1176527719232179e-12
, -3.4345706136796902e-15
, 7.106374596228206e-18
, -9.458532831665803e-21
, 7.309749404022739e-24
, -2.491523307573798e-27
}},{{ -2.8387131384516184e-11
, 1.1390475107328593e-09
, 0.00023120807604118913
, 1.6318549859722082e-07
, -9.119183828082798e-07
, 2.9842650242491955e-06
, -6.791344387053871e-06
, 1.034695467169681e-05
, -1.0424439576423798e-05
, 6.659022387334458e-06
, -2.44126279038688e-06
, 3.91027749688124e-07
},{ 7.521485645570732e-09
, -1.816652149483637e-08
, 0.0002312443867203917
, -9.340867538111298e-09
, -3.900193886716387e-08
, -6.066157447064198e-10
, 1.0054985854188853e-10
, -7.413693279600719e-12
, 4.326075133873619e-13
, -1.6732651151191353e-14
, 3.6965071571311943e-16
, -3.5373863282348784e-18
},{ 0.004258331526114649
, -0.0010306901220669539
, 0.0003232351088351074
, -3.9691954297397845e-06
, 3.9282978039169386e-08
, -2.583694644208844e-10
, 1.157943624139939e-12
, -3.551824771245545e-15
, 7.334381415854274e-18
, -9.741861567831713e-21
, 7.512790974064114e-24
, -2.55521346873582e-27
}},{{ 3.2593427538476505e-11
, -1.328282549332117e-09
, 0.0002348510882002083
, -1.9594703325374155e-07
, 1.0147468054812727e-06
, -3.677390040416635e-06
, 8.468426618118911e-06
, -1.3047917528158028e-05
, 1.3287629185691009e-05
, -8.576064210072395e-06
, 3.1755570509368616e-06
, -5.13581412247746e-07
},{ -6.913129600654379e-09
, 1.6912245873748825e-08
, 0.00023481237366017085
, 8.928494031354307e-09
, -4.676738036833751e-08
, 5.950619603670335e-10
, -6.177282057651006e-11
, 7.46395120016424e-12
, -4.623208354242776e-13
, 1.7303715031942098e-14
, -3.727226684760003e-16
, 3.534912195545466e-18
},{ 0.004049972049150166
, -0.0010169709797367462
, 0.0003276131429734354
, -4.063973464801853e-06
, 4.0247574512454e-08
, -2.645160099833664e-10
, 1.184032140633662e-12
, -3.626710338126292e-15
, 7.477803806840115e-18
, -9.917014546898889e-21
, 7.63577832196253e-24
, -2.5928591378684755e-27
}},{{ -3.0373913207851176e-11
, 1.2289530608513506e-09
, 0.00023414182179081063
, 1.787142682806643e-07
, -1.0029696112485054e-06
, 3.3148571360618806e-06
, -7.596738438415299e-06
, 1.1655169813402317e-05
, -1.1824612796755025e-05
, 7.606194545372856e-06
, -2.807948905638087e-06
, 4.5289223820045987e-07
},{ -6.360449598983834e-09
, 1.5441935806293193e-08
, 0.0002341468083042324
, 8.047386900342006e-09
, -4.601644262928524e-08
, 5.318088992720956e-10
, -5.318645921864516e-11
, 6.638379170663258e-12
, -4.115831340060668e-13
, 1.5350527514238935e-14
, -3.2970967841361127e-16
, 3.121553382676638e-18
},{ 0.004131277711231571
, -0.0010241422818715498
, 0.00032696217598735737
, -4.047407099423747e-06
, 4.012990541014883e-08
, -2.6417647146910183e-10
, 1.184654315413842e-12
, -3.635413821935335e-15
, 7.509986539326319e-18
, -9.978744150340359e-21
, 7.698068819713321e-24
, -2.6190501528783712e-27
}}};
const double COEFF_VEC_ZD[5][3][12] =  {{{ 4.946362680694488e-09
, -5.39542892680585e-06
, 4425.636011023649
, -0.2681941960942372
, 3417.298363166291
, -431.5192162293103
, 1.0
, 1.0
, 1.0
, 1.0
, 1.0
, 1.0
},{ 4.1011903115451215e-06
, -0.0005115471050440077
, 4425.658452476316
, -0.5832614324402623
, 3407.8283606517807
, -62.39802751881218
, -953.1351049472997
, -610.366097446475
, 116.92203749559378
, 1990.7546475774666
, -1373.075909678475
, -1936.9307648675178
},{ 0.2286523303942496
, -195.7015272701131
, 6373.578667154379
, -7773.504136765312
, 19454.420069572705
, -18401.764600113114
, 10282.774015492796
, -3695.996277846693
, 865.4241891896546
, -127.86391567886841
, 10.838713987985196
, -0.4021133656097916
}},{{ 2.819484521189182e-09
, -4.1251457135546455e-06
, 4423.676811887457
, -0.30012666900355317
, 3421.5721325767668
, -586.115970171017
, 1.0
, 1.0
, 1.0
, 1.0
, 1.0
, 1.0
},{ 3.909164721024655e-06
, -0.00048546807794012627
, 4423.698093802654
, -0.5564692999624266
, 3407.3310620175203
, -62.02072248780034
, -937.3200304899449
, -747.7163695478727
, 685.0688018517071
, 850.2264354985037
, -525.4192273976721
, -1795.5749346734428
},{ -1.6726531006194094
, -176.89314677634115
, 6296.330221063298
, -7610.84521084428
, 19239.378437023744
, -18215.774520592888
, 10175.903458598299
, -3655.045255428165
, 855.1002865956167
, -126.22162202944556
, 10.689396337071239
, -0.39620216575641576
}},{{ 7.821152317675076e-10
, -1.2481232726293996e-06
, 4324.761924877385
, -0.09802945270595535
, 3402.510665957487
, -216.46709591300726
, 1.0
, 1.0
, 1.0
, 1.0
, 1.0
, 1.0
},{ 6.165611723770981e-06
, -0.0007342065939745629
, 4324.79439859197
, -0.7504297550234726
, 3404.74906638286
, -73.28526194832062
, -824.9370998383254
, -754.8921407806326
, 520.3882734435285
, 1015.8798802409349
, -660.707199487619
, -1445.1950916734972
},{ 19.060404414785882
, -406.6294194612978
, 7228.390714228183
, -10047.161363469375
, 22529.61137539918
, -20820.657391700912
, 11498.859240182204
, -4100.505829152638
, 954.3390291002127
, -140.30793659878455
, 11.844475233111996
, -0.4378731615396445
}},{{ 3.575753122801213e-09
, -4.271169812603923e-06
, 4258.417177303658
, -0.2441964656207018
, 3406.038178701222
, -432.116548728145
, 1.0
, 1.0
, 1.0
, 1.0
, 1.0
, 1.0
},{ 4.996978385070158e-06
, -0.0006262472627192571
, 4258.44557043322
, -0.7287860430562283
, 3399.394951383096
, -82.16382133881086
, -700.0471325224597
, -1009.0621917166847
, 888.2075864461431
, 916.5629527603509
, -428.01371765595553
, -2671.8538059519096
},{ 34.79440636150863
, -578.1959273972336
, 7918.547161513377
, -11800.918499212285
, 24846.62414248433
, -22616.831795217448
, 12392.311031267956
, -4395.179185640819
, 1018.6486295490508
, -149.2523993631088
, 12.563387068209334
, -0.46330570935891546
}},{{ 5.42693785711986e-09
, -6.376990142719727e-06
, 4270.550168719809
, -0.3559110563188649
, 3406.8796783435064
, -599.9067022759518
, 1.0
, 1.0
, 1.0
, 1.0
, 1.0
, 1.0
},{ 5.703899871260042e-06
, -0.0007174540358857511
, 4270.582357920091
, -0.8446985582690526
, 3394.706316716718
, -96.30175122062587
, -648.0555562390314
, -1176.0831334780562
, 976.4247111959703
, 1303.0655398309416
, -689.0568208233668
, -3322.358160543462
},{ 29.221079775392464
, -516.9741588288789
, 7658.368719921819
, -11162.990245236448
, 23984.123598388567
, -21943.351194978703
, 12056.379247256837
, -4284.274655381598
, 994.4423641309751
, -145.88704614249554
, 12.293090003766364
, -0.4537523034160383
}}};

#endif