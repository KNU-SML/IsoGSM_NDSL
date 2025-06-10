      subroutine wrdlen
 
      common /hrdwrd/ nbytw,nbitw,nrev,iord(8)
      common /charac/ iascii,iatoe(0:255),ietoa(0:255)
      common /quiet / iprt
 
      character*8 cint,dint
      equivalence (cint,int)
      equivalence (dint,jnt)
      logical     print
 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
 
      print = nbytw.eq.0 .and. iprt.eq.1
 
c  count the bits in a word - max 64 allowed
c  -----------------------------------------
 
      int = 1
      do i=1,65
      int = ishft(int,1)
      if(int.eq.0) goto 10
      enddo
10    if(i.ge.65)       goto 900
      if(mod(i,8).ne.0) goto 901
      nbitw = i
      nbytw = i/8
 
c  index the byte storage order -  high byte to low byte
c  -----------------------------------------------------
 
      jnt = 0
      do i=1,nbytw
      int = ishft(1,(nbytw-i)*8)
      do j=1,nbytw
      if(cint(j:j).ne.dint(j:j)) goto 20
      enddo
20    if(j.gt.nbytw) goto 902
      iord(i) = j
      enddo
 
c  set the noreverse flag - 0=noreverse;1=reverse
c  ----------------------------------------------
 
      nrev = 0
      do i=1,nbytw
      if(iord(i).ne.i) nrev = 1
      enddo
 
c  setup an ascii/ebcdic transaltor and determine which is native
c  --------------------------------------------------------------
 
 
      if(iupm('A',8).eq. 65) then
         iascii = 1
      elseif(iupm('A',8).eq.193) then
         iascii = 0
      else
         call abt('wrdlen - cant determine native language')
      endif
 
      do i=0,255
      ietoa(i) = 0
      iatoe(i) = 0
      enddo
 
      ietoa(  1) =   1
      iatoe(  1) =   1
      ietoa(  2) =   2
      iatoe(  2) =   2
      ietoa(  3) =   3
      iatoe(  3) =   3
      ietoa(  5) =   9
      iatoe(  9) =   5
      ietoa(  7) = 127
      iatoe(127) =   7
      ietoa( 11) =  11
      iatoe( 11) =  11
      ietoa( 12) =  12
      iatoe( 12) =  12
      ietoa( 13) =  13
      iatoe( 13) =  13
      ietoa( 14) =  14
      iatoe( 14) =  14
      ietoa( 15) =  15
      iatoe( 15) =  15
      ietoa( 16) =  16
      iatoe( 16) =  16
      ietoa( 17) =  17
      iatoe( 17) =  17
      ietoa( 18) =  18
      iatoe( 18) =  18
      ietoa( 19) =  19
      iatoe( 19) =  19
      ietoa( 22) =   8
      iatoe(  8) =  22
      ietoa( 24) =  24
      iatoe( 24) =  24
      ietoa( 25) =  25
      iatoe( 25) =  25
      ietoa( 29) =  29
      iatoe( 29) =  29
      ietoa( 31) =  31
      iatoe( 31) =  31
      ietoa( 34) =  28
      iatoe( 28) =  34
      ietoa( 37) =  10
      iatoe( 10) =  37
      ietoa( 38) =  23
      iatoe( 23) =  38
      ietoa( 39) =  27
      iatoe( 27) =  39
      ietoa( 45) =   5
      iatoe(  5) =  45
      ietoa( 46) =   6
      iatoe(  6) =  46
      ietoa( 47) =   7
      iatoe(  7) =  47
      ietoa( 50) =  22
      iatoe( 22) =  50
      ietoa( 53) =  30
      iatoe( 30) =  53
      ietoa( 55) =   4
      iatoe(  4) =  55
      ietoa( 60) =  20
      iatoe( 20) =  60
      ietoa( 61) =  21
      iatoe( 21) =  61
      ietoa( 63) =  26
      iatoe( 26) =  63
      ietoa( 64) =  32
      iatoe( 32) =  64
      ietoa( 74) =  91
      iatoe( 91) =  74
      ietoa( 75) =  46
      iatoe( 46) =  75
      ietoa( 76) =  60
      iatoe( 60) =  76
      ietoa( 77) =  40
      iatoe( 40) =  77
      ietoa( 78) =  43
      iatoe( 43) =  78
      ietoa( 79) =  33
      iatoe( 33) =  79
      ietoa( 80) =  38
      iatoe( 38) =  80
      ietoa( 90) =  93
      iatoe( 93) =  90
      ietoa( 91) =  36
      iatoe( 36) =  91
      ietoa( 92) =  42
      iatoe( 42) =  92
      ietoa( 93) =  41
      iatoe( 41) =  93
      ietoa( 94) =  59
      iatoe( 59) =  94
      ietoa( 95) =  94
      iatoe( 94) =  95
      ietoa( 96) =  45
      iatoe( 45) =  96
      ietoa( 97) =  47
      iatoe( 47) =  97
      ietoa(106) = 124
      iatoe(124) = 106
      ietoa(107) =  44
      iatoe( 44) = 107
      ietoa(108) =  37
      iatoe( 37) = 108
      ietoa(109) =  95
      iatoe( 95) = 109
      ietoa(110) =  62
      iatoe( 62) = 110
      ietoa(111) =  63
      iatoe( 63) = 111
      ietoa(121) =  96
      iatoe( 96) = 121
      ietoa(122) =  58
      iatoe( 58) = 122
      ietoa(123) =  35
      iatoe( 35) = 123
      ietoa(124) =  64
      iatoe( 64) = 124
      ietoa(125) =  39
      iatoe( 39) = 125
      ietoa(126) =  61
      iatoe( 61) = 126
      ietoa(127) =  34
      iatoe( 34) = 127
      ietoa(129) =  97
      iatoe( 97) = 129
      ietoa(130) =  98
      iatoe( 98) = 130
      ietoa(131) =  99
      iatoe( 99) = 131
      ietoa(132) = 100
      iatoe(100) = 132
      ietoa(133) = 101
      iatoe(101) = 133
      ietoa(134) = 102
      iatoe(102) = 134
      ietoa(135) = 103
      iatoe(103) = 135
      ietoa(136) = 104
      iatoe(104) = 136
      ietoa(137) = 105
      iatoe(105) = 137
      ietoa(145) = 106
      iatoe(106) = 145
      ietoa(146) = 107
      iatoe(107) = 146
      ietoa(147) = 108
      iatoe(108) = 147
      ietoa(148) = 109
      iatoe(109) = 148
      ietoa(149) = 110
      iatoe(110) = 149
      ietoa(150) = 111
      iatoe(111) = 150
      ietoa(151) = 112
      iatoe(112) = 151
      ietoa(152) = 113
      iatoe(113) = 152
      ietoa(153) = 114
      iatoe(114) = 153
      ietoa(161) = 126
      iatoe(126) = 161
      ietoa(162) = 115
      iatoe(115) = 162
      ietoa(163) = 116
      iatoe(116) = 163
      ietoa(164) = 117
      iatoe(117) = 164
      ietoa(165) = 118
      iatoe(118) = 165
      ietoa(166) = 119
      iatoe(119) = 166
      ietoa(167) = 120
      iatoe(120) = 167
      ietoa(168) = 121
      iatoe(121) = 168
      ietoa(169) = 122
      iatoe(122) = 169
      ietoa(173) =  91
      iatoe( 91) = 173
      ietoa(176) =  48
      iatoe( 48) = 176
      ietoa(177) =  49
      iatoe( 49) = 177
      ietoa(178) =  50
      iatoe( 50) = 178
      ietoa(179) =  51
      iatoe( 51) = 179
      ietoa(180) =  52
      iatoe( 52) = 180
      ietoa(181) =  53
      iatoe( 53) = 181
      ietoa(182) =  54
      iatoe( 54) = 182
      ietoa(183) =  55
      iatoe( 55) = 183
      ietoa(184) =  56
      iatoe( 56) = 184
      ietoa(185) =  57
      iatoe( 57) = 185
      ietoa(189) =  93
      iatoe( 93) = 189
      ietoa(192) = 123
      iatoe(123) = 192
      ietoa(193) =  65
      iatoe( 65) = 193
      ietoa(194) =  66
      iatoe( 66) = 194
      ietoa(195) =  67
      iatoe( 67) = 195
      ietoa(196) =  68
      iatoe( 68) = 196
      ietoa(197) =  69
      iatoe( 69) = 197
      ietoa(198) =  70
      iatoe( 70) = 198
      ietoa(199) =  71
      iatoe( 71) = 199
      ietoa(200) =  72
      iatoe( 72) = 200
      ietoa(201) =  73
      iatoe( 73) = 201
      ietoa(208) = 125
      iatoe(125) = 208
      ietoa(209) =  74
      iatoe( 74) = 209
      ietoa(210) =  75
      iatoe( 75) = 210
      ietoa(211) =  76
      iatoe( 76) = 211
      ietoa(212) =  77
      iatoe( 77) = 212
      ietoa(213) =  78
      iatoe( 78) = 213
      ietoa(214) =  79
      iatoe( 79) = 214
      ietoa(215) =  80
      iatoe( 80) = 215
      ietoa(216) =  81
      iatoe( 81) = 216
      ietoa(217) =  82
      iatoe( 82) = 217
      ietoa(224) =  92
      iatoe( 92) = 224
      ietoa(226) =  83
      iatoe( 83) = 226
      ietoa(227) =  84
      iatoe( 84) = 227
      ietoa(228) =  85
      iatoe( 85) = 228
      ietoa(229) =  86
      iatoe( 86) = 229
      ietoa(230) =  87
      iatoe( 87) = 230
      ietoa(231) =  88
      iatoe( 88) = 231
      ietoa(232) =  89
      iatoe( 89) = 232
      ietoa(233) =  90
      iatoe( 90) = 233
      ietoa(240) =  48
      iatoe( 48) = 240
      ietoa(241) =  49
      iatoe( 49) = 241
      ietoa(242) =  50
      iatoe( 50) = 242
      ietoa(243) =  51
      iatoe( 51) = 243
      ietoa(244) =  52
      iatoe( 52) = 244
      ietoa(245) =  53
      iatoe( 53) = 245
      ietoa(246) =  54
      iatoe( 54) = 246
      ietoa(247) =  55
      iatoe( 55) = 247
      ietoa(248) =  56
      iatoe( 56) = 248
      ietoa(249) =  57
      iatoe( 57) = 249
 
c  show some results
c  -----------------
 
      if(print) then
         print100,nbytw,nbitw,nrev,(iord(i),i=1,nbytw)
         if(iascii.eq.0) print*,'ebcdic is native language'
         if(iascii.eq.1) print*,'ascii  is native language'
      endif
100   format(' wrdlen:nbytw=',i1,' nbitw=',i2,' irev=',i1,' iord=',8i1)
 
      return
900   call abt('wrdlen - a word is more than 64 bits')
901   call abt('wrdlen - a word is not made of bytes')
902   call abt('wrdlen - byte order checking mistake')
      end
