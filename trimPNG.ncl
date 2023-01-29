;=======================================================================
;=======================================================================
undef("trimPNG")
procedure trimPNG(fig_file)
local res,cmd
begin
    strver = systemfunc("ncl -V") 
    ver = tofloat(strver)
    sub = str_sub_str(strver,ver+".","")
    if .not.isStrSubset(fig_file,".png") then 
        fig_file = fig_file + ".png" 
    end if
    cmd = "convert -trim "+fig_file+" "+fig_file
    if (ver.ge.6.2) .or. ((ver.eq.6.2).and.(sub.ge.1)) then
        if fileexists(fig_file) then system(cmd) end if
    else
        if isfilepresent(fig_file) then system(cmd) end if
    end if
    print("")
    print("    "+fig_file)
    print("")
end
;=======================================================================
;=======================================================================