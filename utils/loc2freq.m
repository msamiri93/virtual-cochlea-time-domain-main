function freq = loc2freq(x)
% %  function freq = loc2freq(x)
% %  freq in [kHz]
% %  loc in [mm]
% %  Gerbil freq-location relation from Greenwwod, 1991
% %
% %  YJ Liu

xum = x*1000;
freq = 0.4*(10.^((12e3-xum)*0.174*1e-3)-0.85);
