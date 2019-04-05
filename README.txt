This file describes how to use the Matlab aircraft sizing and wing design scripts for Part 2 of the UEFC Project. 

Aircraft sizing
---------------

The sizing optimization is an extension of the capability developed in Part 1.  These updated scripts include better physical modeling (e.g. Reynolds number and thickness dependence of profile drag; dependence of max thrust on velocity; impact of circular motion on span efficiency; better accounting for dihedral; etc).  While you can certainly look at any of the provided Matlab functions, the functions you are most likely to use are briefly described below.

GetUEFC: As before, this function sets most of the parameters of the design space.  A new parameter from the previous Part 1 versions is CLmax. By setting this, you will constrain the optimization to make sure that CL < CLmax.  You need to do this in order to ensure that the sectional cl stays below the stall limit.  See the slides from the discussion in class for more about that.

GetWpay: Sets the payload weight.  The current version of the script already includes the burrito mass.  However, it does not included any packaging, mounts, etc you may use.  So, you will probably need to update this function as you determine how you will carry the payload.

GetCDpay: Sets the payload-dependent drag coefficient increment.  Currently, this is set to zero (i.e. there is no drag caused by the payload).  Clearly, this is almost certainly incorrect and you can include a payload drag increment here.

scan_V: Essentially the same scan_V as before.  This is the code that will search over (AR,S) to determine optimal designs for each AR,S considered.  Then, it plots contours of V, N, CL, T, d/b.  As well, this new version of scan_V place an asterick (*) at the location in (AR,S) which has the highest velocity.  Finally, scan_V also now prints out what the performance, operating conditions, weight breakdowns, etc are for this highest velocity aircraft.

opt_V(AR,S): this function determines the maximum flight speed achievable for an airplane with the inputted values of AR, S.  This function is called repeatedly by scan_V as it scans over (AR,S).

report_opt_V(AR,S): this function is a wrapper for opt_V.  Calling it will printout the optimized performance, operating conditions, etc found after running opt_V.  It calls opt_V for you and then prints out useful information.

Wing design
-----------

You need to use this script to design the wing twist to achieve good performance (i.e. high span efficiency) while minimizing the potential for bad stall behavior.  The specific pros & cons of wing performance are discussed in the detail notes available on the Stellar site.

The only function you will directly interact with for performing the detailed wing design is UEFC_wvl.  'wvl' stands for Weissinger Vortex Lattice method, which is the underlying aerodynamic model being used to estimate the performance of the wing.

Here are the inputs and outputs to that script:

[e0, maxcl] = UEFC_wvl(AR,S,agr,agt,CL,iplot)

Note: UEFC_wvl call GetUEFC to get lambda, dihedral.  In addition, when you call UEFC_wvl, you supply it with the following inputs:

S: wing area
AR: wing aspect ratio
agr: geometric twist angle at root (deg)
agt: geometric twist angle at tip (deg)
CL: desired wing CL
iplot: plots lift and cl distribution if true

IMPORTANT: agr and agt are really the only parameters you control to improve the aerodynamic performance including stall behavior for a given AR,S,lambda,dihedral, and CL.

Outputs from UEFC_wvl:

e0: span efficiency (then CDi = CL^2/(pi*AR*e0) assuming straight flight
maxcl: maximum sectional lift coefficient

Note that e0 output from UEFC_wvl should be then used to get e0 in GetUEFC.

Also, if you find your maxcl is larger than required, and that you are unable to achieve a maxcl that is below the requirement to avoid stall, then you will need to change your sizing optimization, most likely by lowering the CLmax value until the wing is at a CL for which you can achieve sectional cl's that are below your limit.

As well, if iplot = 1 (true), then plots of lift (ccl) and sectional cl (cl) as well as the planform geometry will be displayed.


