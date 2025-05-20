*Reference: Martinez, L. R. (2022). How Much Should We Trust the Dictator's GDP Growth Estimates? Journal of Political Economy, 130(10).
*Link to the article: <https://doi.org/10.1086/720458>
*https://www.journals.uchicago.edu/doi/suppl/10.1086/720458/suppl_file/20190733data.zip

*Replication file prepared by: Ali Kaan AKSIT, Merve YILMAZ, Tristan GOURMELEN R.

*Academic year 2024/2025

*Replicating : Figure 1, Table 1 column 3, and Figure A8

* The following packages are required to run this do-file
ssc install outreg2
ssc install estout
ssc install binscatter
ssc install parmest
ssc install kountry,replace
ssc install splitvallabels,replace
ssc install reghdfe,replace
ssc install xtabond2,replace
ssc install winsor2,replace
ssc install spmap,replace
ssc install geo2xy
ssc install outreg2
ssc install wbopendata
ssc install labutil
net install http://fmwww.bc.edu/RePEc/bocode/i/ineqdec0


* Before running this do-file, replace the global "dir" in line 19 and the "cd" directory in line 27.

clear all
set more off
global dir "XXX"	

global tables "${dir}\Output\Tables"
global figures "${dir}\Output\Figures"
global temp "${dir}\Data\temp"
global sigmas "${dir}\Data\temp\sigmas"


cd "XXX"

use ZAM_sample_in.dta, clear

gen ctryf = countrycode
encode ctryf, gen(ctry_f) 
gen year_f = year

regress lngdp14 lndn13 fiw fiw2 lndn13_fiw i.ctry_f i.year_f, vce(robust)

estimates store reg_model
outreg2 using regression_results.doc, replace word // Requires `outreg2` package for exporting.


regress lngdp14 lndn13 fiw fiw2 lndn13_fiw i.ctry_f i.year

predict cooksd, cooksd    // Cook's Distance
predict dfits, dfits      // DFITS

gen threshold_cooksd = 4 / e(N) // Common threshold for Cook's Distance
gen out_zam_sign = (cooksd > threshold_cooksd)

gen threshold_dfits = 2 * sqrt(e(k) / e(N)) // DFITS threshold
gen out_zam_sig = (abs(dfits) > threshold_dfits)

gen out_zam_both = (out_zam_sign == 1 & out_zam_sig == 1)

regress lngdp14 lndn13 fiw fiw2 lndn13_fiw i.ctry_f i.year, vce(robust)

save ZAM_sample_out.dta, replace




global seed 20051702	//seed used for bootstrap in estimation of sigma SE
global reps=1000	//number of repetitions for bootstrap

local i=1	// this local identifies each estimate of sigma

timer clear 1
timer on 1

// Figure 1
use "${dir}\Data\master\Estimations", clear

preserve
replace gdp14_growth=gdp14_growth*100
replace dn13_growth=dn13_growth*100
sum dn13_growth if year>=1992&year<=2013&fiw!=.,det
binscatter gdp14_growth dn13_growth if dn13_growth>=-30&dn13_growth<=50&year>=1992&year<=2013&fiw!=.,by(autocracyFH) colors(gs10 gs5) msymbols(circle_hollow triangle) linetype(lfit) legend(r(1) order(2 "Autocracy" 1 "Democracy") pos(11) ring(0)) ytitle("GDP growth (%)") ylabel(,angle(horizontal)) xtitle("Growth of Lights Digital Number (%)") xlabel(-30(10)40) scheme(sj) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) name(raw_scatter,replace) reportreg
graph export "${figures}\figure_1.pdf", as(eps) replace
restore

* *




* *

// Table 1

local out "outreg2 using "${tables}\table_1_c3.tex", bracket nocons nor2 keep(lndn* fiw* Dpfree Dnfree autocracyFH) dec(3)"

preserve

*column 3
capture program drop sigma_boot
program define sigma_boot, rclass
	reghdfe lngdp14 lndn13 fiw lndn13_fiw, absorb(countrycode year) cluster(countrycode)
	sum fiw if e(sample),det
	local iqr: di %2.1f r(p75)-r(p25)
	return scalar sigma = `iqr'*_b[lndn13_fiw]/_b[lndn13]
exit
end
bootstrap r(sigma), cluster(countrycode) idcluster(ctrycode) reps(${reps}) seed($seed): sigma_boot
estat bootstrap,bc p
mat a = e(b)
local sigma: di %4.3f a[1,1]
mat b = e(V)
local sigmaSE: di %4.3f sqrt(b[1,1])
mat c = e(ci_bc)
local sigma_95_low: di %3.2f c[1,1]
local sigma_95_high: di %3.2f c[2,1]
parmest,saving("${sigmas}\sigma`i'",replace) idnum(`i') idstr(T1C3) erows(ci_percentile ci_bc)
reghdfe lngdp14 lndn13 fiw lndn13_fiw, absorb(countrycode year) cluster(countrycode)
`out' addtext(Country FE, Yes, Year FE, Yes, sigma,`sigma', sigma SE, [`sigmaSE'], sigma 95 CI, [`sigma_95_low' `sigma_95_high']) addstat(Countries, e(N_clust), R-squared, e(r2_within))
local i=`i'+1

restore




// Figure A8, robustness test



*create rounded version of FIW, create dummies and interact w/ lights
capture confirm variable fiw_round
if _rc != 0 { 
  
    gen fiw_round = round(fiw)
} 
else {
   
    replace fiw_round = round(fiw)
}

forvalues x=0/6 {
    capture confirm variable Dfiw_`x'
    if _rc != 0 {
        gen Dfiw_`x' = (fiw_round == `x') if fiw_unadj != .
    } 
    else {
        replace Dfiw_`x' = (fiw_round == `x') if fiw_unadj != .
    }

    capture confirm variable lndn13_Dfiw_`x'
    if _rc != 0 {
        gen lndn13_Dfiw_`x' = lndn13 * Dfiw_`x'
    } 
    else {
        replace lndn13_Dfiw_`x' = lndn13 * Dfiw_`x'
    }
}

*run regression for figure
reghdfe lngdp14 lndn13_Dfiw_0 lndn13_Dfiw_1 lndn13_Dfiw_2 lndn13_Dfiw_3 lndn13_Dfiw_4 lndn13_Dfiw_5 lndn13_Dfiw_6,absorb(countrycode i.year fiw_round) cluster(countrycode) nocons
parmest,saving("${temp}\fiw_levels",replace)
*bootstrapped estimates of sigma for bottom
capture program drop sigma_boot
program define sigma_boot, rclass
	reghdfe lngdp14 lndn13_Dfiw_0 lndn13_Dfiw_1 lndn13_Dfiw_2 lndn13_Dfiw_3 lndn13_Dfiw_4 lndn13_Dfiw_5 lndn13_Dfiw_6,absorb(countrycode i.year fiw_round) cluster(countrycode) nocons
	forvalues x=1/6{
		return scalar sigma`x' = (_b[lndn13_Dfiw_`x']-_b[lndn13_Dfiw_0])/_b[lndn13_Dfiw_0]
	}
exit
end
bootstrap r(sigma1) r(sigma2) r(sigma3) r(sigma4) r(sigma5) r(sigma6), cluster(countrycode) idcluster(ctrycode) reps(${reps}) seed($seed): sigma_boot
estat bootstrap,bc p
mat a = e(b)
forvalues x=1/6{
	local sigma`x': di %4.3f a[1,`x']	
}
mat b = e(V)
forvalues x=1/6{
	local sigmaSE`x': di %4.3f sqrt(b[`x',`x'])	
}
*draw graph
preserve
use "${temp}\fiw_levels",clear
gen fiw=_n-1
twoway (bar estimate fiw if fiw==0, barwidth(0.8) color(green%70)) (bar estimate fiw if fiw==1, barwidth(0.8) color(green%30)) (bar estimate fiw if fiw==2, barwidth(0.8) color(yellow%20)) (bar estimate fiw if fiw==3, barwidth(0.8) color(yellow%50)) (bar estimate fiw if fiw==4, barwidth(0.8) color(orange%80)) (bar estimate fiw if fiw==5, barwidth(0.8) color(orange_red%80)) (bar estimate fiw if fiw==6, barwidth(0.8) color(red%90)) (rcap min95 max95 fiw ,lcolor(gs8)), ytitle(NTL Elasticity of GDP) ylabel(0(0.1)0.5, angle(horizontal)) xtitle(Freedom in the World Index) xlabel(0 1 `" `"1"' `"{&sigma}: `sigma1'"'  `"(`sigmaSE1')"' "' 2 `" `"2"' `"{&sigma}: `sigma2'"'  `"(`sigmaSE2')"' "' 3 `" `"3"' `"{&sigma}: `sigma3'"'  `"(`sigmaSE3')"' "' 4 `" `"4"' `"{&sigma}: `sigma4'"'  `"(`sigmaSE4')"' "' 5 `" `"5"' `"{&sigma}: `sigma5'"'  `"(`sigmaSE5')"' "' 6 `" `"6"' `"{&sigma}: `sigma6'"'  `"(`sigmaSE6')"' "' , labsize(small)) legend(off) scheme(sj) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) plotregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white)) name(fh_levels,replace)
graph export "${figures}\figure_A8.pdf", as(pdf) replace
restore
