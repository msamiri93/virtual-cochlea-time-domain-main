function loc = freq2loc(freq)
% %  function loc = freq2loc(freq)
% %  freq in [kHz]
% %  loc in [mm]
% %  Gerbil freq-location relation from Greenwwod, 1991
% %
% %  YJ Liu

xum = 12e3 -(log10(freq/0.4+0.85))/0.174e-3;

loc = xum*1e-3;