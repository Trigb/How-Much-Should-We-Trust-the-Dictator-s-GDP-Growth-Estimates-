
*ssc install outreg2
*ssc install estout
*ssc install binscatter

clear all
set more off
global dir "C:\Users\trist\OneDrive\Documents\Cours\Master Economie du Developpement\M2 EDD\Econometrics Seminar\Replication\rawdata"	

global tables "${dir}\Output\Tables"
global figures "${dir}\Output\Figures"
global temp "${dir}\Data\temp"
global sigmas "${dir}\Data\temp\sigmas"

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
graph export "${figures}\figure_1.eps", as(eps) replace
restore

* *




* *

// Table 1
use "${dir}\Data\master\Estimations", clear

local out "outreg2 using "${tables}\table_1_c3.tex", bracket nocons nor2 keep(lndn* fiw* Dpfree Dnfree autocracyFH) dec(3)"

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





// Figure A8, robustness test

use "${dir}\Data\master\Estimations", clear

*create rounded version of FIW, create dummies and interact w/ lights
gen fiw_round=round(fiw)
forvalues x=0/6{
	gen Dfiw_`x'=(fiw_round==`x') if fiw_unadj!=.
	gen lndn13_Dfiw_`x'=lndn13*Dfiw_`x'
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
