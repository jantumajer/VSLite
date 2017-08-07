pkg load all
more off

%%%%% Parametrizace Bayes modelu pro nalezeni parametru T1, T2, M1, M2 pro jednotlive plochy


a=[1,2,3,4]
calstart=1; calend=73

%for i = [21,22,23]
for i = [21,22,23,26,33,35,42,44,47,48,61,66,76,81,100,106,107,118,124,125,126,127,128,138,140,150,165,189,192,213,231,239,259,292,294,301,307,342,344,349,361,365,373,376,383,390,391,425,436,437,444,491,499,501,506,512,534,539,541,548,566,574,600,610,613,616,633,634,640,667,668,675,678,679,703,706,712,729,762,765,796,804,818,819,826,860,887,888,896,904,938,948,990,992,996,997,1062,1100,1121,1134,1158,1163,1180,1203,1204,1211,1214,1232,1263,1288,1310,1333,1334,1356,1377,1381,1420,1429,1448,1475,1482,1498,1532,1555,1569,1576,1577,1584,1587]
%  Nacteni vstupnich promennych pro jednotlive plochy
%  eval(["load chronologie.mat " (["a", num2str(i)])])
  eval(["load chronologie_zscore.mat " (["a", num2str(i)])])
  eval(["load phi.mat " (["a", num2str(i),"phi"])])
  eval(["load Teplota.mat " (["a", num2str(i),"t"])])
  eval(["load Srazky.mat " (["a", num2str(i),"p"])])
  
%  Model - Bayes
  [T1,T2,M1,M2] = estimate_vslite_params_v2_3(eval(["a", num2str(i),"t"])(:,calstart:calend),eval(["a", num2str(i),"p"])(:,calstart:calend),eval(["a", num2str(i),"phi"]),eval(["a", num2str(i)])(calstart:calend,1)','nsamp',2000, 'intwindow', [-4 12])
  disp(i)
  b=[T1,T2,M1,M2]
  a=[a;b]
    clear T1 T2 M1 M2 b
    eval(["clear chronologie.mat " (["a", num2str(i)])])
    eval(["clear phi.mat " (["a", num2str(i),"phi"])])
    eval(["clear Teplota.mat " (["a", num2str(i),"t"])])
    eval(["clear Srazky.mat " (["a", num2str(i),"p"])])
  
  save Bayes.mat a
  endfor
 
 
%%%%% Vypocet modelovanych letokruhovych krivek pro jednotlive plochy - original podle STW

syear=1; eyear=73
c=[,1940:2012]
d=[0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24]'
x=[0;1]'
for i = [21,22,23,26,33,35,42,44,47,48,61,66,76,81,100,106,107,118,124,125,126,127,128,138,140,150,165,189,192,213,231,239,259,292,294,301,307,342,344,349,361,365,373,376,383,390,391,425,436,437,444,491,499,501,506,512,534,539,541,548,566,574,600,610,613,616,633,634,640,667,668,675,678,679,703,706,712,729,762,765,796,804,818,819,826,860,887,888,896,904,938,948,990,992,996,997,1062,1100,1121,1134,1158,1163,1180,1203,1204,1211,1214,1232,1263,1288,1310,1333,1334,1356,1377,1381,1420,1429,1448,1475,1482,1498,1532,1555,1569,1576,1577,1584,1587]

  eval(["load parametry_z1_73.mat " (["a", num2str(i),"t1"])])
  eval(["load parametry_z1_73.mat " (["a", num2str(i),"t2"])])
  eval(["load parametry_z1_73.mat " (["a", num2str(i),"m1"])])
  eval(["load parametry_z1_73.mat " (["a", num2str(i),"m2"])])
  %  pripadne je mozne pouzit i jiny typ datoveho souboru (1_26, 1_50, z1_26 nebo z1_50)
  
  eval(["load Teplota.mat " (["a", num2str(i),"t"])])
  eval(["load Srazky.mat " (["a", num2str(i),"p"])])
  eval(["load phi.mat " (["a", num2str(i),"phi"])])
  
  [trw_model,gT,gM,gE,M] = VSLite_v2_3(syear,eyear,eval(["a", num2str(i),"phi"]),eval(["a", num2str(i),"t1"]),eval(["a", num2str(i),"t2"]),eval(["a", num2str(i),"m1"]),eval(["a", num2str(i),"m2"]),eval(["a", num2str(i),"t"]),eval(["a", num2str(i),"p"]))
  disp(i)
  
  c=[c;trw_model]
    eval(["clear parametry_z1_73.mat " (["a", num2str(i),"t1"])])
    eval(["clear parametry_z1_73.mat " (["a", num2str(i),"t2"])])
    eval(["clear parametry_z1_73.mat " (["a", num2str(i),"m1"])])
    eval(["clear parametry_z1_73.mat " (["a", num2str(i),"m2"])])
    eval(["clear Teplota.mat " (["a", num2str(i),"t"])])
    eval(["clear Srazky.mat " (["a", num2str(i),"p"])])
    eval(["clear phi.mat " (["a", num2str(i),"phi"])])
    clear trw_model M
    
  g=[repmat(eval(num2str(i)), 73, 1)]; j=[repmat(eval(num2str(i)), 12, 1)]
  
  h=[g,gT',gM']
  d=[d;h]
  
  k=[j,gE]
  x=[x;k]
    clear e gE gM gT
  save Model.mat c
  save GrowthResponse.mat d
  save GR_daylength.mat x
endfor


%%%%% Vypocet modelovanych letokruhovych krivek pro jednotlive plochy - uprava podle Acevedo et al. (2016)

syear=1; eyear=73
c=[,1940:2012]
d=[0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24]'
x=[0;1]'
for i = [21,22,23,26,33,35,42,44,47,48,61,66,76,81,100,106,107,118,124,125,126,127,128,138,140,150,165,189,192,213,231,239,259,292,294,301,307,342,344,349,361,365,373,376,383,390,391,425,436,437,444,491,499,501,506,512,534,539,541,548,566,574,600,610,613,616,633,634,640,667,668,675,678,679,703,706,712,729,762,765,796,804,818,819,826,860,887,888,896,904,938,948,990,992,996,997,1062,1100,1121,1134,1158,1163,1180,1203,1204,1211,1214,1232,1263,1288,1310,1333,1334,1356,1377,1381,1420,1429,1448,1475,1482,1498,1532,1555,1569,1576,1577,1584,1587]

  eval(["load parametry_z1_73.mat " (["a", num2str(i),"t1"])])
  eval(["load parametry_z1_73.mat " (["a", num2str(i),"t2"])])
  eval(["load parametry_z1_73.mat " (["a", num2str(i),"m1"])])
  eval(["load parametry_z1_73.mat " (["a", num2str(i),"m2"])])
  %  pripadne je mozne pouzit i jiny typ datoveho souboru (1_26, 1_50, z1_26 nebo z1_50)
  
  eval(["load Teplota.mat " (["a", num2str(i),"t"])])
  eval(["load Srazky.mat " (["a", num2str(i),"p"])])
  eval(["load phi.mat " (["a", num2str(i),"phi"])])
  
  [trw_model,gT,gM,gE,M] = VSLite_v2_3_Acevedo(syear,eyear,eval(["a", num2str(i),"phi"]),eval(["a", num2str(i),"t1"]),eval(["a", num2str(i),"t2"]),eval(["a", num2str(i),"m1"]),eval(["a", num2str(i),"m2"]),eval(["a", num2str(i),"t"]),eval(["a", num2str(i),"p"]))
  disp(i)
  
  c=[c;trw_model]
    eval(["clear parametry_z1_73.mat " (["a", num2str(i),"t1"])])
    eval(["clear parametry_z1_73.mat " (["a", num2str(i),"t2"])])
    eval(["clear parametry_z1_73.mat " (["a", num2str(i),"m1"])])
    eval(["clear parametry_z1_73.mat " (["a", num2str(i),"m2"])])
    eval(["clear Teplota.mat " (["a", num2str(i),"t"])])
    eval(["clear Srazky.mat " (["a", num2str(i),"p"])])
    eval(["clear phi.mat " (["a", num2str(i),"phi"])])
    clear trw_model M
    
  g=[repmat(eval(num2str(i)), 73, 1)]; j=[repmat(eval(num2str(i)), 12, 1)]
  
  h=[g,gT',gM']
  d=[d;h]
  
  k=[j,gE]
  x=[x;k]
    clear e gE gM gT
  save Model.mat c
  save GrowthResponse.mat d
  save GR_daylength.mat x
endfor