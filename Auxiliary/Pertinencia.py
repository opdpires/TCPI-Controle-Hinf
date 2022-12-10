for i in range(0,32):
    str = "{0:05b} - ".format(i)
    bitHi = int(str[0]) +1
    bitHj = int(str[1]) +1
    bitHk = int(str[2]) +1
    bitHp = int(str[3]) +1
    bitHq = int(str[4]) +1
    
    formato = "h({}) = Hi({})*Hj({})*Hk({})*Hp({})*Hq({});".format(i+1,bitHi,bitHj,bitHk,bitHp,bitHq)
    print(formato)
