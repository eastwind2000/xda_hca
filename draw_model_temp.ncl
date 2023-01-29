
; load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl" 
; load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_inv_string.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

load "./trimPNG.ncl"

begin

; ======================================================================
; cdate  =  "2020030520"
; fh_bk           = (/"036", "048", "060", "072" /)
; ======================================================================


fig_dir ="../figs/"

nvars = 1
nlevs = 1

flist = systemfunc("ls -1  /mnt/i/fcst2023/coldwave/g3km/20230114/*20230114000000*.grb2")
; /mnt/i/fcst2023/coldwave/g3km/20230114/Z_NAFP_C_BABJ_20230114000000_P_NWPC-GRAPES-3KM-ORIG-00000.grb2

nt=dimsizes(flist)

do iregion = 1, 1
    
    if (iregion .eq. 0) then
        ;  全国 region00
        regionx = "region00"
        minlat = 20
        maxlat = 45
        minlon = 100
        maxlon = 135
        ; ref_veclength =  0.012
        ; pts_skip = 5
    end if 

    if (iregion .eq. 1) then
        ; focus fine region
        regionx = "region01"
        minlat = 30
        maxlat = 35
        minlon = 115
        maxlon = 123

    end if 

; =================================================================================== Starting Plot

  
    print("====== Step5: Drawing figs:  " + " ======")

    
    plev = 85000

    do k = 0, 24 ;nt-1   ; time-dimension
;============================================================== getting ecmwf file =======================================================

        ; ListSetType(era5_flist, "join")

        fname  = flist(k)

        print(fname)

        txtbufr = tochar(systemfunc("basename " + fname))

        cdateinfo = tostring( txtbufr(9:)  )

        ; print(cdate)
        ; sleep(200)

        fs00     = addfile(fname, "r")
       
        hlev   = fs00->HGT_P0_L100_GLL0(  {plev}, {minlat:maxlat}, {minlon:maxlon} )
        ulev   = fs00->UGRD_P0_L100_GLL0( {plev}, {minlat:maxlat}, {minlon:maxlon} )
        vlev   = fs00->VGRD_P0_L100_GLL0( {plev}, {minlat:maxlat}, {minlon:maxlon} )
        qlev   = fs00->SPFH_P0_L100_GLL0( {plev}, {minlat:maxlat}, {minlon:maxlon} )
        
        tlev  = fs00->TMP_P0_L100_GLL0(  {plev}, {minlat:maxlat}, {minlon:maxlon} ) -273.15
        copy_VarCoords(hlev, tlev)        
        printVarSummary(tlev)

        model_name = "CMA-MESO Temperature "   

        fig_prefix = fig_dir + "CMA_MESO" + cdateinfo      ;+ tostring(k)

        print(fig_prefix)

        print("===========================================")
        
        ;
        ; lat    = fs00->latitude
        ; lon    = fs00->longitude
        ; copy_VarCoords(h500_v012, h500_bias)       
        ; printVarSummary(t2m_bias)
        ; printVarSummary(datafcst)
        ; sleep(200)
        ;


        spd = sqrt(ulev^2.0 + vlev^2.0)
        copy_VarCoords(ulev, spd)

        ; lat = ulev&latitude
        ; lon = vlev&longitude
        ; databufr = uv2dv_cfd(u850(::3, ::3), v850(::3, ::3), u850&latitude(::3), u850&longitude(::3), 0)*1e6        
        ; div850   = smth9(databufr, 0.25, 0.25, False)
        ; div850!0 = "latitude"
        ; div850!1 = "longitude"
        ; div850&latitude = lat(::3)
        ; div850&longitude = lon(::3)
        ; copy_VarCoords(u850, div850)
        ; timeinfo = grbfs[:]->initial_time0
        ; printVarSummary(prcp)

        print("==========================================================================")
        

        print("Generating figname: " + fig_prefix + ".eps" )

        wks_type                                = "png"
        wks_type@wkWidth                        = 3200
        wks_type@wkHeight                       = 2500
        
        ; wks_type@wkOrientation                 = "landscape"

        wks = gsn_open_wks(wks_type, fig_prefix)   

        ; ; ==================================================================== draw map


        mpres                                   = True
        mpres@gsnDraw                           = False
        mpres@gsnFrame                          = False
        mpres@gsnMaximize                       = False
       

        mpres@mpLimitMode                       = "Corners"
        mpres@mpLeftCornerLatF                  = minlat
        mpres@mpLeftCornerLonF                  = minlon
        mpres@mpRightCornerLatF                 = maxlat
        mpres@mpRightCornerLonF                 = maxlon

        mpres@mpOutlineOn                       = False
        mpres@mpFillOn                          = False

        mpres@pmTickMarkDisplayMode             =  "Always"
        mpres@tmXBLabelFontHeightF              = 0.012

        mpres@gsnLeftString                     = model_name + " " + cdateinfo + " " + tostring(plev)        
        mpres@gsnLeftStringOrthogonalPosF       = 0
        mpres@gsnLeftStringFontHeightF          = 0.01

        mpres@gsnCenterString                   = ""

        ; mpres@gsnRightStringFontHeightF         = 0.008

        plot_mp  = gsn_csm_map(wks, mpres)

        getvalues plot_mp                                   ; important for overplay control
        "vpXF"      : xpos
        "vpYF"      : ypos
        "vpWidthF"  : width
        "vpHeightF" : height
        end getvalues

        print("xpos = " + xpos)

        ; ; ==================================================================== Adding shapefiles
        
        shp_path = "/home/lse/res/ncarg/shp/"
        shp_list =(/ "bou1_4l.shp", "continents_lines.shp", "hyd1_4l.shp", "province_l.shp" /)
        shp_color = (/"Black",      "Black",                "Blue",        "Black"  /)
        nshp       = 4
        shpplot    = new(nshp, graphic)
        shpres                  = True
        do ishp = 0, nshp - 1 
            shpres@gsLineColor      = shp_color(ishp)
            shpres@gsLineThicknessF = 1.5
            shpplot(ishp) = gsn_add_shapefile_polylines(wks, plot_mp, shp_path + shp_list(ishp), shpres)  
        end do

        ; ==============================================================================plot spd

        biasres                                   =  True
        biasres@vpXF                              =  xpos
        biasres@vpYF                              =  ypos
        biasres@vpWidthF                          =  width
        biasres@vpHeightF                         =  height

        biasres@gsnDraw                           = False
        biasres@gsnFrame                          = False
        biasres@gsnAddCyclic                      = False
        biasres@gsnMaximize                       = False       ; pay attention to correcness of overlay
        biasres@gsnLeftStringFontHeightF          = 0.01
        biasres@gsnRightStringFontHeightF         = 0.01
        biasres@gsnCenterStringFontHeightF        = 0.01
        biasres@gsnCenterString                   = ""
            
        biasres@cnFillMode                        = "CellFill"

        ; biasres@trGridType                        = "TriangularMesh"
        ; biasres@cnFillMode                        = "RasterFill"        
        ; biasres@cnLevelSelectionMode              = "ExplicitLevels"    ; 
        ; biasres@cnLevels                          = (/ 0.1, 1, 5, 10, 25, 30, 40, 50, 75, 100, 125, 150, 175, 200/)

        biasres@cnLevelSelectionMode              = "ManualLevels"    ;      
        biasres@cnMaxLevelValF                    =   2
        biasres@cnMinLevelValF                    =  -2 
        biasres@cnLevelSpacingF                   =  0.25

        biasres@cnFillOn                          = True
        cmap                                      = read_colormap_file("temp_diff_18lev")
        biasres@cnFillPalette                     = cmap(:, :)
    
        biasres@cnLineLabelInterval               = 1       
        biasres@cnLinesOn                         = False
        biasres@cnLineLabelsOn                    = False  
        biasres@cnLineThicknessF                  = 4.0

        biasres@cnInfoLabelOn                     = False
        biasres@pmLabelBarDisplayMode             = "Always"
        ; spdres@pmLabelBarSide                    = "Right"

        biasres@lbLabelStride                     =  2
        biasres@lbLabelFontHeightF                =  0.0015

        biasres@pmLabelBarOrthogonalPosF          =  0.01       ; Move labelbar up
        biasres@pmLabelBarParallelPosF            =  0.5        ; Move labelbar up   
        biasres@pmLabelBarHeightF                 =  0.4
        biasres@pmLabelBarWidthF                  =  0.04


        biasres@tmXBLabelsOn                      = False
        biasres@tmYLLabelsOn                      = False        

        ; biasres@cnFillDrawOrder                   = "Predraw"
        ; spdres@lbOverrideFillOpacity             = True
        ; spdres@lbFillOpacityF                    = 0
        
        biasres@lbLabelBarOn                      = False   ; importrant settings when overlay map !

        ; ==============================================================================plot hgt
        
        t2mres                                   =  True

        t2mres@vpXF                              =  xpos
        t2mres@vpYF                              =  ypos
        t2mres@vpWidthF                          =  width
        t2mres@vpHeightF                         =  height

        t2mres@gsnDraw                           = False
        t2mres@gsnFrame                          = False
        t2mres@gsnAddCyclic                      = False
        t2mres@gsnMaximize                       = False        ; pay attention to correcness of overlay 

        t2mres@gsnLeftString                     = ""
        t2mres@gsnRightString                    = ""

        ; t2mres@trGridType                        = "TriangularMesh"        

        t2mres@cnLevelSelectionMode              = "ManualLevels"    ; 
        t2mres@cnMaxLevelValF                    =  20
        t2mres@cnMinLevelValF                    = -20
        t2mres@cnLevelSpacingF                   =  1

        t2mres@cnFillOn                          = True
        cmap_temp                                = read_colormap_file("ViBlGrWhYeOrRe")
        t2mres@cnFillPalette                     = cmap_temp
        
        ; t2mres@cnLevelSelectionMode              = "ExplicitLevels"
        ; t2mres@cnLevels                          = (/ 0 /)

        t2mres@cnLineLabelInterval               = 1
        t2mres@cnLineLabelBackgroundColor        = -1
        t2mres@cnLineLabelFontHeightF            = 0.01


        t2mres@cnLinesOn                         = True        
        t2mres@cnLineThicknessF                  = 3.0
        t2mres@cnLineColor                       = "gray80"

        ; t2mres@gsnContourZeroLineThicknessF      = 5.0


        t2mres@cnInfoLabelOn                     = "False"
        
        t2mres@tmXTOn                            = False
        t2mres@tmYLOn                            = False
        t2mres@tmXBLabelsOn                      = False       
        t2mres@tmYLLabelsOn                      = False

        ; t2mres@pmTickMarkDisplayMode             =  "Always"

        t2mres@lbLabelBarOn                      = False

        t2mres@pmLabelBarDisplayMode             = "Always"
        t2mres@pmLabelBarOrthogonalPosF          =  0.04       ; Move labelbar up
        t2mres@pmLabelBarParallelPosF            =  0.5        ; Move labelbar up   
        t2mres@pmLabelBarHeightF                 =  0.4
        t2mres@pmLabelBarWidthF                  =  0.06
        

        ;=================================================================================== 


        tmzres                                   =  True

        tmzres@vpXF                              =  xpos
        tmzres@vpYF                              =  ypos
        tmzres@vpWidthF                          =  width
        tmzres@vpHeightF                         =  height

        tmzres@gsnDraw                           = False
        tmzres@gsnFrame                          = False
        tmzres@gsnAddCyclic                      = False
        tmzres@gsnMaximize                       = False        ; pay attention to correcness of overlay 

        tmzres@gsnLeftString                     = ""
        tmzres@gsnRightString                    = ""

        ; tmzres@trGridType                        = "TriangularMesh"        

        tmzres@cnLevelSelectionMode              = "ExplicitLevels"    ; 
        tmzres@cnLevels                          =  (/0/)

        tmzres@cnFillOn                          = False

        tmzres@cnLineLabelInterval               = 1
        tmzres@cnLineLabelBackgroundColor        = -1
        tmzres@cnLineLabelFontHeightF            = 0.01


        tmzres@cnLinesOn                         = True        
        tmzres@cnLineThicknessF                  = 4.0
        tmzres@cnLineColor                       = "red"

        ; tmzres@gsnContourZeroLineThicknessF      = 5.0


        tmzres@cnInfoLabelOn                     = "False"
        
        tmzres@tmXTOn                            = False
        tmzres@tmYLOn                            = False
        tmzres@tmXBLabelsOn                      = False       
        tmzres@tmYLLabelsOn                      = False

        ; tmzres@pmTickMarkDisplayMode             =  "Always"

        tmzres@lbLabelBarOn                      = False
        

        ;=================================================================================== 

        vecres = True
        
        vecres@vpXF                        = xpos
        vecres@vpYF                        = ypos
        vecres@vpWidthF                    =  width
        vecres@vpHeightF                   = height

        vecres@gsnLeftString               = " "
        vecres@gsnRightString              = " "

        vecres@gsnDraw                    = False
        vecres@gsnFrame                   = False
        vecres@gsnMaximize                = False
        vecres@gsnAddCyclic               = False

        vecres@vcGlyphStyle                = "WindBarb"
        vecres@vcWindBarbScaleFactorF      =  2.5       ; important settings !

        ; vecres@vcWindBarbColor             = "blue4"
        vecres@vcMonoWindBarbColor        = False
        vecres@vcLevelSelectionMode       = "ExplicitLevels"
        vecres@vcLevels                   = (/   12,  16, 18, 20 /)
        vecres@vcLevelColors              = (/"blue4", "brown1", "brown2","brown3", "brown4" /)

        vecres@vcWindBarbLineThicknessF    =  2.0

        vecres@vcWindBarbTickLengthF       =  0.35
        vecres@vcWindBarbTickSpacingF      =  0.14
        vecres@vcWindBarbCalmCircleSizeF   =  0.01

        vecres@vcRefAnnoOn                =  False

        vecres@vcRefMagnitudeF            =  10
        vecres@vcRefLengthF               = 0.015      ; length reference vector
        vecres@vcMinDistanceF             = 0.015

        vecres@vcMinDistanceF             = 0.015

        vecres@tmXTOn                            = False
        vecres@tmYLOn                            = False
        vecres@tmXBLabelsOn                      = False       
        vecres@tmYLLabelsOn                      = False

        ; ==================================================================================
        ; vecres@vcRefAnnoOrthogonalPosF    =  0.05
        ; vecres@vcLineArrowColor           =  "blue4"       
        ; =================================================================================== 


        ; plot_qpe  = gsn_csm_contour(wks, databufr( {minlat:maxlat}, {minlon:maxlon} ), qperes)

        ; plot_bias = gsn_csm_contour(wks, t2m_bias, biasres )

        plot_wind = gsn_csm_vector(wks, ulev, vlev, vecres)
        
        plot_t2m  = gsn_csm_contour(wks, tlev, t2mres )

        plot_tmz  = gsn_csm_contour(wks, tlev, tmzres )
        ; plot_slp  = gsn_csm_contour(wks, slp(k, :, :), slpres )
        
    ; ========================================================   Adding track and marks

        ; plres                              = True
        ; plres@gsLineThicknessF             = 6.0
        ; plres@gsLineColor                  = "blue"
        ; trackplot                          = gsn_add_polyline(wks, plot_mp, cen_lon, cen_lat, plres)

        ; tfDoNDCOverlay 

        overlay(plot_mp, plot_t2m)

        overlay(plot_mp, plot_tmz)   

        overlay(plot_mp, plot_wind)        

        ; overlay(plot_mp, plot_slp)               
     
        ;=================================================================================== 
        
        maximize_output(wks, True)                       ; important function for overlay 

        ; delete(plot_hgt)
        delete(plot_mp)
        ; delete(plot_slp)
        delete(shpplot)
        delete(shpres)
        ; biasres@cnFillOn                          = True
        ; delete(slpres)
        delete(t2mres)
        delete(mpres)
        delete(wks)      

        ; delete(prcp)

        delete(ulev)
        delete(vlev)
        ; delete(lat)
        ; delete(lon)        

        trimPNG(fig_prefix+".png")

    end do ; for k level
    
end do
; ==================================================================== Preforme interpolation using ESMF

end




