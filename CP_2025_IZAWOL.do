*******************************************************************************
** Replication code for Cattaneo and Palomba (2025)
** "Leveraging Covariates in Regression Discontinuity Designs"
** Last modified: 2025-04-21
*******************************************************************************


** Install packages
* net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace
* net install rdhte, from(https://raw.githubusercontent.com/rdpackages/rdhte/main/stata) replace

clear all
set more off

**************************************
* 1) Set up paths and load data
**************************************
use "headstart.dta", clear

**************************************
* 2) Generate analysis variables
**************************************
* outcome, running, heterogeneity covariate
global y mort_age59_related_postHS
global x povrate60
gen byte census1960_pop_ind = (census1960_pop >= 10000)

* cutoff
global c = 59.1968

* covariates for efficiency
global z "census1960_pctblack census1960_pctsch1417 census1960_pctsch534 census1960_pctsch25plus census1960_pop1417 census1960_pop534 census1960_pop25plus census1960_pcturban"

* covariate for heterogeneity
global w_hte "census1960_pop_ind"

**************************************
* 3) Compute optimal bandwidths
**************************************
* all units
quietly: rdrobust $y $x, c($c) kernel(triangular) p(1)
scalar h_all = e(h_l)

* by subgroup w==0 and w==1
forvalues g = 0/1 {
    quietly: rdrobust $y $x if $w_hte==`g', c($c) kernel(triangular) p(1)
    scalar h`g' = e(h_l)
}


***********************************************
* 4) Figure 1: global and local RD (all units)
***********************************************

cap drop rdplot_*
rdplot $y $x, c($c) p(1) kernel("tri") genvars hide

* global
twoway (scatter rdplot_mean_y rdplot_mean_bin if inrange($x,40,70), msize(small) mcolor(black)) ///
(function `e(eq_l)', range(40 $c) lcolor(black) lwidth(thick) lpattern(solid)) 			 		///
(function `e(eq_r)', range($c 70) lcolor(black) lwidth(thick) lpattern(solid)), 		 		///
xline($c, lcolor(black) lwidth(medthin)) legend(off) ytitle("child mortality") xtitle("poverty rate")

* zoomed
global c_l = $c - h_all
global c_r = $c + h_all
twoway (scatter rdplot_mean_y rdplot_mean_bin if inrange($x,40,70), msize(small) mcolor(gray))  ///
(scatter rdplot_mean_y rdplot_mean_bin if inrange($x,$c_l,$c_r), msize(small) mcolor(black)) 	///
(function `e(eq_l)', range($c_l $c) lcolor(black) lwidth(thick) lpattern(solid)) 			 	///
(function `e(eq_r)', range($c $c_r) lcolor(black) lwidth(thick) lpattern(solid)), 		 	 	///
xline($c_l, lcolor(black) lwidth(medthin)) xline($c_r, lcolor(black) lwidth(medthin)) 			///
xline($c, lcolor(black) lwidth(medthin)) legend(off) ytitle("child mortality") xtitle("poverty rate")



********************************************
* 5) Figure 2a & 2b: heterogeneous RD by w
********************************************
cap drop rdplot_*
rdplot $y $x if $w_hte == 0, c($c) p(1) kernel("tri") genvars hide

matrix coef_r0 = e(coef_r)
matrix coef_l0 = e(coef_l)

local eq_l0 = "  y = coef_l0[1, 1]*(x-59.1968)^0   + coef_l0[1+1, 1]*(x-59.1968)^1"
local eq_r0 = "  y = coef_r0[1, 1]*(x-59.1968)^0   + coef_r0[1+1, 1]*(x-59.1968)^1"

foreach var of varlist rdplot_* {
	rename `var' w0_`var'
}

rdplot $y $x if $w_hte == 1, c($c) p(1) kernel("tri") genvars hide

foreach var of varlist rdplot_* {
	rename `var' w1_`var'
}

* global
twoway (scatter w0_rdplot_mean_y w0_rdplot_mean_bin if inrange($x,40,70), msize(small) mcolor(red)) ///
(scatter w1_rdplot_mean_y w1_rdplot_mean_bin if inrange($x,40,70), msize(small) mcolor(green))      ///
(function `eq_l0', range(40 $c) lcolor(red) sort lwidth(thick) lpattern(solid)) 			 		///
(function `eq_r0', range($c 70) lcolor(red) sort lwidth(thick) lpattern(solid)) 		 			///
(function `e(eq_l)', range(40 $c) lcolor(green) sort lwidth(thick) lpattern(solid)) 			 	///
(function `e(eq_r)', range($c 70) lcolor(green) sort lwidth(thick) lpattern(solid)), 		 		///
legend(order(1 2) row(1) pos(6) lab(1 "population <10k") lab(2 "population >10k")) 					///
xline($c, lcolor(black) lwidth(medthin)) ytitle("child mortality") xtitle("poverty rate")


* zoomed
global c_l_0 = $c - h0
global c_r_0 = $c + h0
global c_l_1 = $c - h1
global c_r_1 = $c + h1

twoway (scatter w0_rdplot_mean_y w0_rdplot_mean_bin if inrange($x,40,70), msize(small) mcolor(gray))   ///
(scatter w1_rdplot_mean_y w1_rdplot_mean_bin if inrange($x,40,70), msize(small) mcolor(gray)) 		   ///
(scatter w0_rdplot_mean_y w0_rdplot_mean_bin if inrange($x,$c_l_0,$c_r_0), msize(small) mcolor(red))   ///
(scatter w1_rdplot_mean_y w1_rdplot_mean_bin if inrange($x,$c_l_1,$c_r_1), msize(small) mcolor(green)) ///
(function `eq_l0', range($c_l_0 $c) lcolor(red) sort lwidth(thick) lpattern(solid)) 			 	   ///
(function `eq_r0', range($c $c_r_0) lcolor(red) sort lwidth(thick) lpattern(solid)) 		 		   ///
(function `e(eq_l)', range($c_l_1 $c) lcolor(green) sort lwidth(thick) lpattern(solid)) 			   ///
(function `e(eq_r)', range($c $c_r_1) lcolor(green) sort lwidth(thick) lpattern(solid)), 		 	   ///
legend(order(3 4) row(1) pos(6) lab(3 "population <10k") lab(4 "population >10k")) 					   ///
xline($c_l_0, lcolor(red) lwidth(medthin)) xline($c_r_0, lcolor(red) lwidth(medthin)) 				   ///
xline($c_l_1, lcolor(green) lwidth(medthin)) xline($c_r_1, lcolor(green) lwidth(medthin)) 			   ///
xline($c, lcolor(black) lwidth(medthin)) ytitle("child mortality") xtitle("poverty rate")


********************************************
* 5) Table 1: covs for efficiency
********************************************

* Check "population" is pre-treatment
rdplot $w_hte $x, c($c)
rdrobust $w_hte $x, c($c) rho(1) vce("hc3")

* RD on all units w/o covs
rdrobust $y $x, c($c) rho(1) vce("hc3")

* RD on all units w covs for efficiency
rdrobust $y $x, c($c) rho(1) vce("hc3") covs($z)

* RD on all units w all covs for efficiency
rdrobust $y $x, c($c) rho(1) vce("hc3") covs($z census1960_pop)


********************************************
* 5) Table 2: heterogeneity analysis 
********************************************
** Panel A) different bandwidths
rdhte $y $x, c($c) covs_hte(census1960_pop_ind) vce("hc3")
rdhte_lincom 0.census1960_pop_ind - 1.census1960_pop_ind
rdhte $y $x, c($c) covs_hte(census1960_pop_ind) vce("hc3") covs_eff($z)
rdhte_lincom 0.census1960_pop_ind - 1.census1960_pop_ind

** Panel B) same bandwidths
rdhte $y $x, c($c) covs_hte(census1960_pop_ind) vce("hc3") bwjoint
rdhte_lincom 0.census1960_pop_ind - 1.census1960_pop_ind
rdhte $y $x, c($c) covs_hte(census1960_pop_ind) vce("hc3") covs_eff($z) bwjoint
rdhte_lincom 0.census1960_pop_ind - 1.census1960_pop_ind


























