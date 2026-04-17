program racialgap

use toolbox

implicit none

integer, parameter :: nage = 66
integer, parameter :: nret = 46
integer, parameter :: na = 250 ! asset grid size
integer, parameter :: ns = 3 ! employed v.s. unemployed
integer, parameter :: nr = 2 ! balck and white
integer, parameter :: niab  = 3 ! human capital grid size
integer, parameter :: neta  = 5 ! income shock grid size
! integer, parameter :: mid_eta  = 2 ! income shock grid size
integer, parameter :: maxiter= 1000
integer, parameter :: useold = 1 ! =1: import initial guess

real*8, parameter :: agrow = 0.06d0
real*8, parameter :: penkappa = 0.4d0 ! pension replacement rate
real*8, parameter :: ww = 0.5d0 ! converge speed of K
real*8, parameter :: ww1 = 0.2d0 ! converge speed of crate
real*8, parameter :: amax = 130.0d0
real*8, parameter :: amin = 0.0d0
real*8, parameter :: tolk  = 1d-5
real*8, parameter :: tolc  = 1d-3
real*8, parameter :: beta = 0.989d0 ! time discount
real*8, parameter :: alfa = 0.64d0 ! labor share in production function
real*8, parameter :: dep = 0.08d0 ! depreciation rate
real*8, parameter :: piep = 0.1717d0 ! probabiliry the criminal is caught (pi_a) 1996-2011
real*8, parameter :: sigma1 = 1.5d0 ! risk aversion
real*8, parameter :: bar = 1.0d0    ! TFP
real*8, parameter :: tau1 = 0.911d0, tau2 = 0.053d0 ! tax function parameter
real*8, parameter :: n_p =  0.01d0 ! population growth
real*8, parameter :: nu = 0.24d0 ! share of accidental bequest goes to newborns
real*8, parameter :: kappa  = 0.0575d0 ! maximum transfer can received
! ******************************
! ******************************
real*8, parameter :: cbar  = 0.098d0 ! consumption for incarcerated person (c_bar)
real*8, parameter :: gamma = 0.064d0 ! fraction of stolen income
real*8, parameter :: hbar  = 0.403d0  ! fraction of time spent on working

real*8 :: a_bor(0:na), u_temp(0:1), drugpro(nage,nr), beqtrans(nage,nr)
real*8 :: y(nage,niab,ns,neta,nr), confis(3), psi(nage+1,niab,nr), erate1(ns,nr), eta_wgt(neta,nr)
real*8 :: cbarage(nage), eta(niab,nr,neta), pieta(niab,nr,neta,neta), steal(nage)
real*8, dimension(nr) :: crrace, beq, first_a, frac_debt, procrrace, drugcrrace, beqdis, phi, uifrac
real*8, dimension(nage,niab,nr) :: eff, muraceedu
real*8, dimension(nage,0:na,niab,ns,neta,nr,5,0:1) :: fdis, aplus, trans, cons, util, income
real*8, dimension(nage,0:na,niab,ns,neta,nr,5,0:1,0:1) :: dis, fdis_new
real*8, dimension(nage,0:na,niab,ns,neta,nr) :: allinc, alltaxinc, inctax, ev, pi_v
real*8, dimension(nage*(na+1)*niab*ns*neta*nr) :: sort_inc, dis_relinc, cum_relinc
integer, dimension(nage*(na+1)*niab*ns*neta*nr) :: order_relinc
real*8, dimension(nage,0:na,niab,ns,neta,nr,0:1) :: val, cr
real*8, dimension(nage,ns,ns,nr) :: pis
real*8, dimension(niab,nr):: ylevel, poplevel, pen, rho, sigma_eps ! persistence/variance of income shock
real*8, dimension(nage,niab,ns,neta,nr,5) :: preincome
real*8, dimension(niab,nr,2) :: ie_pro


real*8 :: ll, mpl, kdev, cdev, rbar, atr, dinc, taup, totalpen, stealinc, prorate, drugrate
real*8 :: exdem, proratenew, x1, x2, x3, neg_util, llnew, taxinc
real*8 :: kknew, inves, ybar, cratenew, tot_debt, unben
real*8 :: aprime, kk, crate, iavg, totinc, workpop, losefrac, scale_factor, loseaux
real*8 :: out, varphi, cons_jail, incar, cc, tottrans, avgtaxinc, incar_w, incar_b
real*8 :: cons_com, efftotal, expense, revenue, taxbase, gg, totinctax
real*8 :: incargroup(10), lambda(nage,ns-1,ns-1,nr), num_cr(2), vic(5), inc_cutoff(4)
real*8 :: pie(nage,0:na,niab,ns,neta,nr,5), ie_aux

integer :: i, j, iter, iage, ia, iab, ii, ial, iar, ir, ip, is, ie, ie_p, is_p, ios, ie_scarring, ie_pp
integer :: i_rec, i_rec_new


!! time variables and open mp
integer :: starttime, rate, endtime
real*8 :: seconds
!!!!!*********************************************************************!!!!!!!

!!!!! measure the time need to run the program
call system_clock(starttime, rate)

!!!!!*********************************************************************!!!!!!!

!!!!!*******probability of transfering to lower state when released********!!!!!!!



ie_pro(1,1,:) = (/0.0d0,1.0d0/)
ie_pro(2,1,:) = (/0.05d0,0.95d0/)
ie_pro(3,1,:) = (/0.0d0,1.0d0/)
ie_pro(1,2,:) = (/0.0d0,1.0d0/)
ie_pro(2,2,:) = (/0.5d0,0.5d0/)
ie_pro(3,2,:) = (/0.0d0,1.0d0/)


!!!!!***********************************************!!!!!



!!!!!*******step-wise victimization probability********!!!!!!!
vic(1) = 4.989210143685341
vic(2) = 4.745407775044441
vic(3) = 4.358159750699997
vic(4) = 2.874727174639702
vic(5) = 2.584454417228699

vic = vic/100.0d0

!!!!!***********************************************!!!!!

!!!!!!!!!!!****************** income shock block *****************!!!!!!
rho(1,1) = 0.9582d0
rho(2,1) = 0.9955d0
rho(3,1) = 0.8422d0
rho(1,2) = 0.9999d0
rho(2,2) = 0.9752d0
rho(3,2) = 0.9817d0
sigma_eps(1,1) = 0.0236d0
sigma_eps(2,1) = 0.0132d0
sigma_eps(3,1) = 0.0456d0
sigma_eps(1,2) = 0.0088d0
sigma_eps(2,2) = 0.0286d0
sigma_eps(3,2) = 0.0307d0


do iab = 1,niab
    do ir = 1,nr
        call discretize_AR(rho(iab,ir), 0d0, sigma_eps(iab,ir), eta(iab,ir,:), pieta(iab,ir,:,:))
    end do
end do

eta = exp(eta)

!!!!!*********************************************************************!!!!!!!


!!!!!!!!!!!****************** life uncertainty block *****************!!!!!!
open(unit=65,file='../model_input/surv_b_3edu.txt',status="old",action="read",position="rewind")
open(unit=66,file='../model_input/surv_w_3edu.txt',status="old",action="read",position="rewind") 
    do iage = 1, nage
        read(unit=65, fmt=*) psi(iage,:,1) 
        read(unit=66, fmt=*) psi(iage,:,2)
    end do
close(unit=65)
close(unit=66)

psi(nage+1,:,:) = 0.0d0

!!!!!*********************************************************************!!!!!!!

!!!!!!!!!!!********************Population structure block *********************!!!!!!

open(unit=31, file = '../model_input/pop_morg.txt', status="old", action="read", position="rewind")
do iage = 1,nage
    read(unit=31,fmt=*)muraceedu(iage,:,1)
end do
do iage = 1,nage
    read(unit=31,fmt=*)muraceedu(iage,:,2)
end do
close(unit=31)


muraceedu = muraceedu/ sum(muraceedu)
workpop = sum(muraceedu(:nret,:,:))
!!!!!*********************************************************************!!!!!!!

!!!!!!!!!*********************** labor efficiency block ******************!!!!!!
eff = 0.0d0
open(unit=81,file='../model_input/eff_black_psid.txt',status="old",action="read",position="rewind")
open(unit=82,file='../model_input/eff_white_psid.txt',status="old",action="read",position="rewind")
    do iage = 1, nret
        ! age-human capital efficiency profiles
        read(unit=81, fmt=*) eff(iage,1,1),eff(iage,2,1),eff(iage,3,1) 
        read(unit=82, fmt=*) eff(iage,1,2),eff(iage,2,2),eff(iage,3,2)
    end do
close(unit=81)
close(unit=82)


efftotal = 0.0d0
do iage = 1, nret
    do iab = 1, niab
        do ir = 1, nr
            efftotal = efftotal + eff(iage,iab,ir) *muraceedu(iage,iab,ir)/workpop
        end do
    end do
end do

eff = eff/efftotal

!!!!!*********************************************************************!!!!!!!


!!!!!!!!!!!****************** employemnt transition block *****************!!!!!!
uifrac = (/0.27d0, 0.39d0/) 
phi = (/0.6978d0, 0.7538d0/) ! unemployment benefit replacement rate

open(unit=91,file='../model_input/urate_black.txt',status="old",action="read",position="rewind")
open(unit=92,file='../model_input/urate_white.txt',status="old",action="read",position="rewind")   
do iage = 1, nage
    read(unit=91, fmt=*) lambda(iage,1,1,1),lambda(iage,1,2,1),lambda(iage,2,1,1),lambda(iage,2,2,1) 
    read(unit=92, fmt=*) lambda(iage,1,1,2),lambda(iage,1,2,2),lambda(iage,2,1,2),lambda(iage,2,2,2) 
end do
close(unit=91)
close(unit=92)


do iage = 1, nage
    do ir = 1, nr
        pis(iage,1,1,ir) = lambda(iage,1,1,ir)
        pis(iage,1,2,ir) = lambda(iage,1,2,ir) * uifrac(ir)
        pis(iage,1,3,ir) = lambda(iage,1,2,ir) * (1-uifrac(ir))
        pis(iage,2,1,ir) = lambda(iage,2,1,ir)
        pis(iage,2,2,ir) = lambda(iage,2,2,ir) * uifrac(ir)
        pis(iage,2,3,ir) = lambda(iage,2,2,ir) * (1-uifrac(ir))
        pis(iage,3,1,ir) = lambda(iage,2,1,ir)
        pis(iage,3,2,ir) = lambda(iage,2,2,ir) * uifrac(ir)
        pis(iage,3,3,ir) = lambda(iage,2,2,ir) * (1-uifrac(ir))
    end do
end do

do iage = 1, nage
    do ir = 1, nr
        pis(iage,1,:,ir)= pis(iage,1,:,ir)/sum((pis(iage,1,:,ir)))
        pis(iage,2,:,ir)= pis(iage,2,:,ir)/sum((pis(iage,2,:,ir)))
        pis(iage,3,:,ir)= pis(iage,3,:,ir)/sum((pis(iage,3,:,ir)))
    end do
end do

erate1(1,1) =  1.0d0- 0.2238136d0          ! employment rate at age 1 of black
erate1(2,1) =  0.2238136d0 * uifrac(1)     ! UI at age 1 of black
erate1(3,1) =  0.2238136d0 * (1-uifrac(1))          
erate1(1,2) =  1.0d0- 0.1080198d0          ! employment rate at age 1 of white
erate1(2,2) =  0.1080198d0 * uifrac(2)     ! UI at age 1 of white
erate1(3,2) =  0.1080198d0 * (1-uifrac(2))   


!!!!!*********************************************************************!!!!!!!

!!!!!!!!!*********************** Exogenous drug crime block ******************!!!!!!
drugpro = 0.0d0
open(unit=22,file='../model_input/drugpro_avg.txt',status="old",action="read",position="rewind")
    do iage = 1, nret
        read(unit=22,fmt=*)drugpro(iage,1)
    end do
    do iage = 1, nret
        read(unit=22,fmt=*)drugpro(iage,2)
    end do
close(unit=22)


!!!!!*********************************************************************!!!!!!!



!!!!!!!!!!!****************** Initial Guess block *****************!!!!!!
if (useold == 1) then
    open(unit=57,file='../model_input/initial_noscarring_drug.txt',status="old",action="read",position="rewind")
        read(unit=57, fmt=*)kk    
        read(unit=57, fmt=*)ll
        read(unit=57, fmt=*)crate 
        read(unit=57, fmt=*)iavg  
        read(unit=57, fmt=*)taup
        read(unit=57, fmt=*)avgtaxinc
        read(unit=57, fmt=*)losefrac
        read(unit=57, fmt=*)ybar
        read(unit=57, fmt=*)prorate
        read(unit=57,fmt=*)beq(1)
        read(unit=57,fmt=*)beq(2)
        do iab = 1, niab
            do ir = 1, nr
                read(unit=57,fmt=*)pen(iab,ir)
            end do
        end do
        do iage = 1, nage
            do ia = 0, na
                do iab = 1, niab
                    do is = 1, ns
                        do ie = 1, neta
                            do ir = 1, nr
                                read(unit=57,fmt=*, iostat=ios)pi_v(iage,ia,iab,is,ie,ir)
                                if (ios /= 0) exit
                            end do
                        end do
                    end do
                end do
            end do
        end do
    close(unit=57)
else
    kk    = 2.07587d0 ! Initial guess of K
    ll    = 0.36095d0
    crate = 0.0408d0
    iavg  = 2.54357d0
    taup  = 0.07570d0
    avgtaxinc = 0.51601d0 
    losefrac = 0.05787080
    ybar = 0.43814d0
    prorate = 0.0352d0
    beqdis = 0.01d0
    beqtrans = 0.1d0
    pi_v = 0.03d0
end if
!!!!!*********************************************************************!!!!!!!

! asset grid
call grid_Cons_Grow(a_bor, amin, amax, agrow)

iter = 1


cbarage = 1d-10
cbarage(1:nret) = cbar



open(unit=10, file="../results/output_noscarring_drug.txt", status="replace", action="write", position="rewind")
write(unit=10,fmt="(4(A20, F8.5))")'gamma= ', gamma, 'cbar= ', cbar, 'alfa= ', alfa
write(unit=10,fmt="(A20,F8.4,A20,I8)")'amax= ', amax, 'na= ', na
write(unit=10,fmt="(4(A20, F8.5))")'agrow= ', agrow, 'dep= ', dep
write(unit=10,fmt="(4(A20, F8.5))")'beta= ', beta,'TFP= ', bar,'pi_a= ', piep
write(unit=10,fmt="(4(A20, F8.5))")'rho (b<hs)= ', rho(1,1), 'rho (b=hs)= ', rho(2,1), 'rho (bco)= ', rho(3,1)
write(unit=10,fmt="(4(A20, F8.5))")'rho (w<hs)= ', rho(1,2), 'rho (w=hs)= ', rho(2,2), 'rho (wco)= ', rho(3,2)
write(unit=10,fmt="(4(A20, F8.5))")'sigma (b<hs)= ', sigma_eps(1,1),'sigma (b=hs)= ', sigma_eps(2,1), 'sigma (bco)= ', sigma_eps(3,1)
write(unit=10,fmt="(4(A20, F8.5))")'sigma (w<hs)= ', sigma_eps(1,2),'sigma (w=hs)= ', sigma_eps(2,2), 'sigma (wco)= ', sigma_eps(3,2)
write(unit=10,fmt="(4(A20, F8.5))")'Frac. work time= ', hbar
write(unit=10,fmt="(4(A20, F8.5))")'Inc Tax Level= ', tau1, 'Inc Tax Progre= ', tau2
write(unit=10,fmt="(4(A20, F8.5))")'SS repl. rate= ', penkappa, 'sigma= ', sigma1
close(unit=10)


do while ( (((kdev > tolk) .or. (cdev > tolc)) .and. (iter < maxiter) ) .or. (iter == 1))

    
    beqtrans = 0.0d0
    beqtrans(31:35,1) = (1.0d0-nu)*beq(1)/sum(muraceedu(31:35,:,1))
    beqtrans(31:35,2) = (1.0d0-nu)*beq(2)/sum(muraceedu(31:35,:,2))

    beqdis(1) = nu*(beq(1))/sum(muraceedu(1,:,1)*(1+n_p))
    beqdis(2) = nu*(beq(2))/sum(muraceedu(1,:,2)*(1+n_p))
    

    pie(:,:,:,:,:,:,1) = 1.0d0-pi_v                 ! 1-pi_v
    pie(:,:,:,:,:,:,2) = pi_v                       ! pi_v
    pie(:,:,:,:,:,:,3) = (1.0d0-piep)*(1.0d0-pi_v)  ! (1-pi_a)*(1-pi_v)
    pie(:,:,:,:,:,:,4) = piep                       ! pi_a
    pie(:,:,:,:,:,:,5) = (1.0d0-piep)*pi_v          ! (1-pi_a)*pi_v

    mpl  = (alfa)*bar*(kk/ll)**(1.0d0-alfa) ! marginal production of labor
    rbar = (1.0d0-alfa)*bar*(kk/ll)**(-alfa) - dep ! interest rate
    atr = 1.0d0 + rbar
    
    steal = gamma * ybar
    steal(nret+1:) = 0.0d0

    cons = 0.0d0
    util = 0.0d0
    aplus = 0.0d0
    val = 0.0d0

    ! calculate age-J agent's policy funtion
    do ia = 0, na
        do iab = 1, niab ! human capital types
            do is = 1, ns
                do ie = 1, neta
                    do ir = 1, nr ! blacks and whites

                        dinc = pen(iab,ir)
                        taxinc = dinc + a_bor(ia)*rbar
                        call incometax(taxinc, x3)
                        inctax(nage,ia,iab,is,ie,ir) = x3
                        totinc = dinc + a_bor(ia)*atr - inctax(nage,ia,iab,is,ie,ir)

                        do ip = 0, 1 ! incarcerated by drug crime
                            if (ip == 0) then ! ip = 0, not incarcerated by drug crime

                                do ii = 1, 5
                                    if (ii == 4) then
                                        trans(nage,ia,iab,is,ie,ir,ii,ip) = 0.0d0
                                    else
                                        call transfer(x1, totinc)
                                        trans(nage,ia,iab,is,ie,ir,ii,ip) = x1
                                    end if
                                end do

                                
                                cons(nage,ia,iab,is,ie,ir,1,ip) = dinc+atr*a_bor(ia)-inctax(nage,ia,iab,is,ie,ir)

                                cons(nage,ia,iab,is,ie,ir,2,ip) = (1-losefrac)*dinc + atr*a_bor(ia) - inctax(nage,ia,iab,is,ie,ir)

                                cons(nage,ia,iab,is,ie,ir,3,ip) = dinc + atr*a_bor(ia) + steal(nage) - inctax(nage,ia,iab,is,ie,ir)

                                cons(nage,ia,iab,is,ie,ir,4,ip) = cbarage(nage)

                                cons(nage,ia,iab,is,ie,ir,5,ip) = (1-losefrac)*dinc + atr*a_bor(ia) + steal(nage) - inctax(nage,ia,iab,is,ie,ir)

                                cons(nage,ia,iab,is,ie,ir,:,ip) = cons(nage,ia,iab,is,ie,ir,:,ip) + trans(nage,ia,iab,is,ie,ir,:,ip)
                                util(nage,ia,iab,is,ie,ir,:,ip) = cons(nage,ia,iab,is,ie,ir,:,ip)**(1-sigma1)/(1-sigma1)

                                aplus(nage,ia,iab,is,ie,ir,:,ip) = 0 
                                
                                u_temp(0) = pie(nage,ia,iab,is,ie,ir,1)*util(nage,ia,iab,is,ie,ir,1,ip) + pie(nage,ia,iab,is,ie,ir,2)*util(nage,ia,iab,is,ie,ir,2,ip) 
                                ! expected utility if choose legal
                                u_temp(1) = pie(nage,ia,iab,is,ie,ir,3)*util(nage,ia,iab,is,ie,ir,3,ip) + pie(nage,ia,iab,is,ie,ir,4)*util(nage,ia,iab,is,ie,ir,4,ip) &
                                + pie(nage,ia,iab,is,ie,ir,5)*util(nage,ia,iab,is,ie,ir,5,ip)
                                ! expeceted utility if choose crime
                                
                                if(u_temp(0) >= u_temp(1)) then
                                    cr(nage,ia,iab,is,ie,ir,ip)  = 0 ! choose legal
                                    val(nage,ia,iab,is,ie,ir,ip) = u_temp(0)
                                else
                                    cr(nage,ia,iab,is,ie,ir,ip)   = 1 ! choose crime
                                    val(nage,ia,iab,is,ie,ir,ip)  = u_temp(1)
                                end if
                                

                            else ! ip = 1, incarcerated by drug crime

                                trans(nage,ia,iab,is,ie,ir,4,ip) = 0.0d0
                                cons(nage,ia,iab,is,ie,ir,4,ip) = cbarage(nage)

                                aplus(nage,ia,iab,is,ie,ir,4,ip) = 0.0d0
                                cr(nage,ia,iab,is,ie,ir,ip)  = 1 ! commit into crime because of drug
                                val(nage,ia,iab,is,ie,ir,ip) = cbarage(nage)**(1-sigma1)/(1-sigma1)

                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do

    call interpolate(nage)

    ! calculate age-1 to age (J-1) agent's policy funtion
    do iage = nage-1, 1, -1
        do ia = 0, na
            do iab = 1, niab ! human capital types
                do is = 1, ns
                    do ie = 1, neta
                        do ir = 1, nr

                            if (iage > nret) then
                                dinc = pen(iab,ir)
                            else
                                if (is == 1) then ! if employed
                                    dinc = mpl * hbar * (1.0d0-taup)*eff(iage,iab,ir)*eta(iab,ir,ie)
                                elseif (is == 2) then ! if unemployed
                                    dinc = mpl * hbar * phi(ir)     *eff(iage,iab,ir)*eta(iab,ir,ie)
                                else
                                    dinc = 0.0d0
                                end if
                            endif
                            taxinc = dinc + a_bor(ia)*rbar
                            call incometax(taxinc, x3)
                            inctax(iage,ia,iab,is,ie,ir) = x3
                            totinc = dinc + a_bor(ia)*atr - inctax(iage,ia,iab,is,ie,ir)+ beqtrans(iage,ir)

                            do ip = 0, 1 ! incarcerated by drug crime
                                if (ip == 0) then ! ip = 0, not incarcerated by drug crime

                                    do ii = 1,5
                                        if (ii == 4) then
                                            trans(iage,ia,iab,is,ie,ir,ii,ip) = 0.0d0
                                            cons(iage,ia,iab,is,ie,ir,ii,ip) = cbarage(iage)
                                            aplus(iage,ia,iab,is,ie,ir,ii,ip) = atr*(a_bor(ia))+ beqtrans(iage,ir)

                                            call linint_Grow(aplus(iage,ia,iab,is,ie,ir,ii,ip), amin, amax, agrow, na, ial, iar, varphi)
                                            ial = max(min(ial, na-1), 0)
                                            iar = max(min(iar, na), 1)
                                            varphi = max(min(varphi, 1.0d0), 0.0d0)
                                            
                                            ! ie_scarring = max(1, ie - 1)
                                            ! util(iage,ia,iab,is,ie,ir,ii,ip) = cbarage(iage)**(1-sigma1)/(1-sigma1) &
                                            ! + beta*psi(iage+1,iab,ir)*(varphi*ev(iage+1,ial,iab,is,ie,ir)+(1-varphi)*ev(iage+1,iar,iab,is,ie,ir))**(1-sigma1)/(1-sigma1)
                                            ! + beta*psi(iage+1,iab,ir)*(varphi*ev(iage+1,ial,iab,is,1,ir)+(1-varphi)*ev(iage+1,iar,iab,is,1,ir))**(1-sigma1)/(1-sigma1)
                                            ! + beta*psi(iage+1,iab,ir)*(varphi*ev(iage+1,ial,iab,is,2,ir)+(1-varphi)*ev(iage+1,iar,iab,is,2,ir))**(1-sigma1)/(1-sigma1)
                                            ! + beta*psi(iage+1,iab,ir)*(varphi*ev(iage+1,ial,iab,is,ie_scarring,ir)+(1-varphi)*ev(iage+1,iar,iab,is,ie_scarring,ir))**(1-sigma1)/(1-sigma1)
                                            
                                            ! + beta*psi(iage+1,iab,ir)*(ie_pro(iab,ir,1)*(varphi*ev(iage+1,ial,iab,is,ie_scarring,ir)+(1-varphi)*ev(iage+1,iar,iab,is,ie_scarring,ir)) &
                                            ! + ie_pro(iab,ir,2)*(varphi*ev(iage+1,ial,iab,is,ie,ir)+(1-varphi)*ev(iage+1,iar,iab,is,ie,ir)))**(1-sigma1)/(1-sigma1)
                                            
                                            ie_aux = 0.0d0
                                            do ie_pp = 1,2
                                                if (ie_pp == 1) ie_scarring = max(1, ie - 1)
                                                if (ie_pp == 2) ie_scarring = ie
                                                ie_aux = ie_aux + ie_pro(iab,ir,ie_pp)*(varphi*ev(iage+1,ial,iab,is,ie_scarring,ir)+(1-varphi)*ev(iage+1,iar,iab,is,ie_scarring,ir))
                                            end do
                                            
                                            util(iage,ia,iab,is,ie,ir,ii,ip) = cbarage(iage)**(1-sigma1)/(1-sigma1) &
                                                + beta*psi(iage+1,iab,ir)*ie_aux**(1-sigma1)/(1-sigma1)
                                        else
                                            call transfer(x1, totinc)
                                            trans(iage,ia,iab,is,ie,ir,ii,ip) = x1

                                            ! get initial guess for the individual choices
                                            x2 = 0.0d0! for very small assets level, search from 0!
                                            
                                            call fminsearch(x2, neg_util, amin, amax, saving)
                                            
                                            aplus(iage,ia,iab,is,ie,ir,ii,ip) = x2
                                            cons(iage,ia,iab,is,ie,ir,ii,ip) = cons_com

                                            call linint_Grow(x2, amin, amax, agrow, na, ial, iar, varphi)
                                            ial = max(min(ial, na-1), 0)
                                            iar = max(min(iar, na), 1)
                                            varphi = max(min(varphi, 1.0d0), 0.0d0)

                                            util(iage,ia,iab,is,ie,ir,ii,ip) = cons(iage,ia,iab,is,ie,ir,ii,ip)**(1-sigma1)/(1-sigma1) &
                                            + beta*psi(iage+1,iab,ir)*(varphi*ev(iage+1,ial,iab,is,ie,ir)+(1-varphi)*ev(iage+1,iar,iab,is,ie,ir))**(1-sigma1)/(1-sigma1)
                                        end if
                                    end do
                                    u_temp(0) = pie(iage,ia,iab,is,ie,ir,1)*util(iage,ia,iab,is,ie,ir,1,ip) + pie(iage,ia,iab,is,ie,ir,2)*util(iage,ia,iab,is,ie,ir,2,ip) 
                                    ! expected utility if choose legal
                                    u_temp(1) = pie(iage,ia,iab,is,ie,ir,3)*util(iage,ia,iab,is,ie,ir,3,ip) + pie(iage,ia,iab,is,ie,ir,4)*util(iage,ia,iab,is,ie,ir,4,ip)&
                                                +pie(iage,ia,iab,is,ie,ir,5)*util(iage,ia,iab,is,ie,ir,5,ip)
                                    ! expeceted utility if choose crime
                                    if(u_temp(0) >= u_temp(1)) then
                                        cr(iage,ia,iab,is,ie,ir,ip)  = 0 ! choose legal
                                        val(iage,ia,iab,is,ie,ir,ip) = u_temp(0)
                                    else
                                        cr(iage,ia,iab,is,ie,ir,ip)   = 1 ! choose crime
                                        val(iage,ia,iab,is,ie,ir,ip)  = u_temp(1)
                                    end if
                                else ! ip = 1, incarcerated by drug crime
                                    trans(iage,ia,iab,is,ie,ir,4,ip) = 0.0d0
                                    cons(iage,ia,iab,is,ie,ir,4,ip) = cbarage(iage) 
                                    aplus(iage,ia,iab,is,ie,ir,4,ip) = atr*a_bor(ia)+ beqtrans(iage,ir)

                                    call linint_Grow(aplus(iage,ia,iab,is,ie,ir,4,ip), amin, amax, agrow, na, ial, iar, varphi)
                                    ial = max(min(ial, na-1), 0)
                                    iar = max(min(iar, na), 1)

                                    cr(iage,ia,iab,is,ie,ir,ip)  = 1 
                                    
                                    ! exp: get rid of scarring for drug
                                    val(iage,ia,iab,is,ie,ir,ip) = cbarage(iage)**(1-sigma1)/(1-sigma1) &
                                    + beta*psi(iage+1,iab,ir)*(ie_pro(iab,ir,1)*(varphi*ev(iage+1,ial,iab,is,ie,ir)+(1-varphi)*ev(iage+1,iar,iab,is,ie,ir)) &
                                    + ie_pro(iab,ir,2)*(varphi*ev(iage+1,ial,iab,is,ie,ir)+(1-varphi)*ev(iage+1,iar,iab,is,ie,ir)))**(1-sigma1)/(1-sigma1)

                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call interpolate(iage)
    end do

    eta_wgt = 0.0d0
    eta_wgt(3,1) = 1.0d0
    eta_wgt(3,2) = 1.0d0


    dis = 0.0d0
    do iab = 1, niab
        do is = 1, ns
            do ie = 1, neta
                do ir = 1, nr
                    call linint_Grow(beqdis(ir) , amin, amax, agrow, na, ial, iar, varphi)
                    ial = max(min(ial, na-1), 0)
                    iar = max(min(iar, na), 1)
                    varphi = max(min(varphi, 1d0), 0d0)

                    dis(1,ial,iab,is,ie,ir,4,1,0) = dis(1,ial,iab,is,ie,ir,4,1,0) + eta_wgt(ie,ir)*drugpro(1,ir)*erate1(is,ir)*varphi

                    dis(1,iar,iab,is,ie,ir,4,1,0) = dis(1,iar,iab,is,ie,ir,4,1,0) + eta_wgt(ie,ir)*drugpro(1,ir)*erate1(is,ir)*(1.0d0-varphi)

                    if (cr(1,ial,iab,is,ie,ir,0) == 0) then
                        do j = 1,2
                            dis(1,ial,iab,is,ie,ir,j,0,0) = dis(1,ial,iab,is,ie,ir,j,0,0) + eta_wgt(ie,ir)*(1-drugpro(1,ir))&
                            *erate1(is,ir)*pie(1,ial,iab,is,ie,ir,j)*varphi
                        end do
                    else
                        do j = 3,5
                            dis(1,ial,iab,is,ie,ir,j,0,0) = dis(1,ial,iab,is,ie,ir,j,0,0) + eta_wgt(ie,ir)*(1-drugpro(1,ir))&
                            *erate1(is,ir)*pie(1,ial,iab,is,ie,ir,j)*varphi
                        end do
                    end if

                    if (cr(1,iar,iab,is,ie,ir,0) == 0) then
                        do j = 1,2
                            dis(1,iar,iab,is,ie,ir,j,0,0) = dis(1,iar,iab,is,ie,ir,j,0,0) + eta_wgt(ie,ir)*(1-drugpro(1,ir))&
                            *erate1(is,ir)*pie(1,iar,iab,is,ie,ir,j)*(1.0d0-varphi)
                        end do
                    else
                        do j = 3,5
                            dis(1,iar,iab,is,ie,ir,j,0,0) = dis(1,iar,iab,is,ie,ir,j,0,0) + eta_wgt(ie,ir)*(1-drugpro(1,ir))&
                            *erate1(is,ir)*pie(1,iar,iab,is,ie,ir,j)*(1.0d0-varphi)
                        end do
                    end if
                end do
            end do
        end do
    end do
    
    do iage = 2, nage
        do ia = 0, na
            do iab = 1, niab
                do is = 1, ns
                    do ie = 1, neta
                        do ir = 1, nr
                            do ii = 1, 5
                                do ip = 0, 1
                                    call linint_Grow(aplus(iage-1,ia,iab,is,ie,ir,ii,ip), amin, amax, agrow, na, ial, iar, varphi)
                                    ial = max(min(ial, na-1), 0)
                                    iar = max(min(iar, na), 1)
                                    varphi = max(min(varphi, 1d0), 0d0)
                                    
                                    do i_rec = 0, 1
                                        
                                        i_rec_new = 0
                                        ! exp: get rid of scarring for drug crime
                                        ! if (ii == 4) then
                                        if (ii == 4 .and. ip == 0) then
                                            i_rec_new = 1
                                            do ie_pp = 1 , 2

                                                if (ie_pp == 1) ie_scarring = max(1, ie - 1)
                                                if (ie_pp == 2) ie_scarring = ie
                                                do is_p = 1, ns
                                                    do ie_p = 1, neta
                                                        dis(iage,ial,iab,is_p,ie_p,ir,4,1,i_rec_new) = dis(iage,ial,iab,is_p,ie_p,ir,4,1,i_rec_new) + drugpro(iage,ir)*&
                                                        pis(iage,is,is_p,ir) *ie_pro(iab,ir,ie_pp)*pieta(iab,ir,ie_scarring,ie_p) * varphi    * dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)
                                                        
                                                        dis(iage,iar,iab,is_p,ie_p,ir,4,1,i_rec_new) = dis(iage,iar,iab,is_p,ie_p,ir,4,1,i_rec_new) + drugpro(iage,ir)*&
                                                        pis(iage,is,is_p,ir) *ie_pro(iab,ir,ie_pp)*pieta(iab,ir,ie_scarring,ie_p) * (1-varphi)* dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)

                                                        
                                                        if (cr(iage,ial,iab,is_p,ie_p,ir,0)==0) then
                                                            do j = 1,2
                                                                dis(iage,ial,iab,is_p,ie_p,ir,j,0,i_rec_new) = dis(iage,ial,iab,is_p,ie_p,ir,j,0,i_rec_new) + pis(iage,is,is_p,ir) &
                                                                *ie_pro(iab,ir,ie_pp)*pieta(iab,ir,ie_scarring,ie_p) *(1-drugpro(iage,ir))*pie(iage,ial,iab,is_p,ie_p,ir,j) * varphi* dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)
                                                            end do
                                                        else
                                                            do j = 3,5
                                                                dis(iage,ial,iab,is_p,ie_p,ir,j,0,i_rec_new) = dis(iage,ial,iab,is_p,ie_p,ir,j,0,i_rec_new) + pis(iage,is,is_p,ir) &
                                                                *ie_pro(iab,ir,ie_pp)*pieta(iab,ir,ie_scarring,ie_p) *(1-drugpro(iage,ir))*pie(iage,ial,iab,is_p,ie_p,ir,j) * varphi* dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)
                                                            end do
                                                        end if

                                                        if (cr(iage,iar,iab,is_p,ie_p,ir,0)==0) then
                                                            do j = 1,2
                                                                dis(iage,iar,iab,is_p,ie_p,ir,j,0,i_rec_new) = dis(iage,iar,iab,is_p,ie_p,ir,j,0,i_rec_new) + pis(iage,is,is_p,ir) &
                                                                *ie_pro(iab,ir,ie_pp)*pieta(iab,ir,ie_scarring,ie_p) *(1-drugpro(iage,ir))*pie(iage,iar,iab,is_p,ie_p,ir,j) * (1-varphi)* dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)
                                                            end do
                                                        else
                                                            do j = 3,5
                                                                dis(iage,iar,iab,is_p,ie_p,ir,j,0,i_rec_new) = dis(iage,iar,iab,is_p,ie_p,ir,j,0,i_rec_new) + pis(iage,is,is_p,ir) &
                                                                *ie_pro(iab,ir,ie_pp)*pieta(iab,ir,ie_scarring,ie_p) *(1-drugpro(iage,ir))*pie(iage,iar,iab,is_p,ie_p,ir,j) * (1-varphi)* dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)
                                                            end do
                                                        end if
                                                    end do
                                                end do
                                            end do
                                        else
                                            ie_scarring = ie;
                                            if (i_rec == 1) then
                                                i_rec_new = 1;
                                            end if

                                            do is_p = 1, ns
                                                do ie_p = 1, neta
                                                    dis(iage,ial,iab,is_p,ie_p,ir,4,1,i_rec_new) = dis(iage,ial,iab,is_p,ie_p,ir,4,1,i_rec_new) + drugpro(iage,ir)*&
                                                    pis(iage,is,is_p,ir) *pieta(iab,ir,ie_scarring,ie_p) * varphi    * dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)
                                                    
                                                    dis(iage,iar,iab,is_p,ie_p,ir,4,1,i_rec_new) = dis(iage,iar,iab,is_p,ie_p,ir,4,1,i_rec_new) + drugpro(iage,ir)*&
                                                    pis(iage,is,is_p,ir) *pieta(iab,ir,ie_scarring,ie_p) * (1-varphi)* dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)

                                                    
                                                    if (cr(iage,ial,iab,is_p,ie_p,ir,0)==0) then
                                                        do j = 1,2
                                                            dis(iage,ial,iab,is_p,ie_p,ir,j,0,i_rec_new) = dis(iage,ial,iab,is_p,ie_p,ir,j,0,i_rec_new) + pis(iage,is,is_p,ir) &
                                                            *pieta(iab,ir,ie_scarring,ie_p) *(1-drugpro(iage,ir))*pie(iage,ial,iab,is_p,ie_p,ir,j) * varphi* dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)
                                                        end do
                                                    else
                                                        do j = 3,5
                                                            dis(iage,ial,iab,is_p,ie_p,ir,j,0,i_rec_new) = dis(iage,ial,iab,is_p,ie_p,ir,j,0,i_rec_new) + pis(iage,is,is_p,ir) &
                                                            *pieta(iab,ir,ie_scarring,ie_p) *(1-drugpro(iage,ir))*pie(iage,ial,iab,is_p,ie_p,ir,j) * varphi* dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)
                                                        end do
                                                    end if

                                                    if (cr(iage,iar,iab,is_p,ie_p,ir,0)==0) then
                                                        do j = 1,2
                                                            dis(iage,iar,iab,is_p,ie_p,ir,j,0,i_rec_new) = dis(iage,iar,iab,is_p,ie_p,ir,j,0,i_rec_new) + pis(iage,is,is_p,ir) &
                                                            *pieta(iab,ir,ie_scarring,ie_p) *(1-drugpro(iage,ir))*pie(iage,iar,iab,is_p,ie_p,ir,j) * (1-varphi)* dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)
                                                        end do
                                                    else
                                                        do j = 3,5
                                                            dis(iage,iar,iab,is_p,ie_p,ir,j,0,i_rec_new) = dis(iage,iar,iab,is_p,ie_p,ir,j,0,i_rec_new) + pis(iage,is,is_p,ir) &
                                                            *pieta(iab,ir,ie_scarring,ie_p) *(1-drugpro(iage,ir))*pie(iage,iar,iab,is_p,ie_p,ir,j) * (1-varphi)* dis(iage-1,ia,iab,is,ie,ir,ii,ip,i_rec)
                                                        end do
                                                    end if
                                                end do
                                            end do
                                        end if

                                        
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    fdis_new = 0.0d0
    do iage = 1, nage
        do iab = 1, niab
            do ir = 1, nr
                fdis_new(iage,:,iab,:,:,ir,:,:,:) = dis(iage,:,iab,:,:,ir,:,:,:)*muraceedu(iage,iab,ir)
            end do
        end do
    end do
    fdis_new = fdis_new/sum(fdis_new)
    if (iter == 1) then
        fdis = sum(fdis_new,dim = 9)
    else
        fdis = (1-ww1)*fdis + ww1*sum(fdis_new,dim = 9)
    end if
    
    y = 0.0d0
    allinc = 0.0d0
    alltaxinc = 0.0d0
    do iage = 1, nage
        do iab = 1, niab
            do ie = 1, neta
                do ir = 1, nr
                    if (iage<=nret)then
                        y(iage,iab,1,ie,ir) = mpl * hbar * (1.0d0-taup) *eff(iage,iab,ir)*eta(iab,ir,ie)
                        y(iage,iab,2,ie,ir) = mpl * hbar * phi(ir)      *eff(iage,iab,ir)*eta(iab,ir,ie)
                        y(iage,iab,3,ie,ir) = 0.0d0
                    else
                        y(iage,iab,:,ie,ir) = pen(iab,ir)
                    endif
                    do ia = 0, na
                        do is = 1, ns
                            alltaxinc(iage,ia,iab,is,ie,ir) = y(iage,iab,is,ie,ir) + a_bor(ia)*rbar
                            allinc(iage,ia,iab,is,ie,ir)    = y(iage,iab,is,ie,ir) + a_bor(ia)*atr - inctax(iage,ia,iab,is,ie,ir)+ beqtrans(iage,ir)
                        end do
                    end do
                end do
            end do
        end do
    end do

    confis = 0.0d0
    cons_jail = 0.0d0

    do iage = 1, nage
        do ia = 0, na
            do iab = 1, niab
                do is = 1, ns
                    do ie = 1, neta
                        do ir = 1, nr
                            confis(3) = confis(3) + fdis(iage,ia,iab,is,ie,ir,4,0) * steal(iage)
                            do ip = 0, 1
                                confis(1) = confis(1) + fdis(iage,ia,iab,is,ie,ir,4,ip) * (y(iage,iab,is,ie,ir) &
                                        - inctax(iage,ia,iab,is,ie,ir)) * (1.0d0 - pi_v(iage,ia,iab,is,ie,ir))
                                confis(2) = confis(2) + fdis(iage,ia,iab,is,ie,ir,4,ip) * (y(iage,iab,is,ie,ir) * (1-losefrac)&
                                        - inctax(iage,ia,iab,is,ie,ir)) * pi_v(iage,ia,iab,is,ie,ir)
                                cons_jail = cons_jail + fdis(iage,ia,iab,is,ie,ir,4,ip) * cbarage(iage)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    incar = sum(fdis(:,:,:,:,:,:,4,:))
    incargroup(1) = sum(fdis(1:5,:,:,:,:,:,4,:))
    do i = 2, 9
        incargroup(i) = sum(fdis((i-1)*5+1:i*5,:,:,:,:,:,4,:))
    end do
    incargroup(10) = sum(fdis(46:nage,:,:,:,:,:,4,:))
    
    incargroup = incargroup*100 / incar 

    incar_b = sum(fdis(:,:,:,:,:,1,4,:))
    incar_w = sum(fdis(:,:,:,:,:,2,4,:))
    do iab = 1, niab
        do ir = 1, nr
            poplevel(iab,ir) = sum(fdis(1:nret,:,iab,:,:,ir,:,:))
        end do
    end do

    first_a = 0.0d0
    proratenew = 0.0d0
    drugrate = 0.0d0
    procrrace = 0.0d0
    drugcrrace = 0.0d0
    cc = 0.0d0
    kknew = 0.0d0
    totinctax = 0.0d0
    tottrans = 0.0d0
    cratenew = 0.0d0
    iavg = 0.0d0
    avgtaxinc = 0.0d0
    ybar = 0.0d0
    num_cr = 0.0d0
    crrace = 0.0d0
    aprime = 0.0d0
    beq = 0.0d0
    ylevel = 0.0d0
    totalpen = 0.0d0
    frac_debt = 0.0d0
    tot_debt = 0.0d0
    do iage = 1, nage
        do ia = 0, na
            do iab = 1, niab
                do is = 1, ns
                    do ie = 1, neta
                        do ir = 1, nr
                            do ii = 1, 5
                                proratenew  = proratenew + cr(iage,ia,iab,is,ie,ir,0) * fdis(iage,ia,iab,is,ie,ir,ii,0) 
                                drugrate    = drugrate   + cr(iage,ia,iab,is,ie,ir,1) * fdis(iage,ia,iab,is,ie,ir,ii,1)
                                procrrace(ir)  = procrrace(ir)  + cr(iage,ia,iab,is,ie,ir,0)  * fdis(iage,ia,iab,is,ie,ir,ii,0) /sum(muraceedu(:,:,ir))
                                drugcrrace(ir) = drugcrrace(ir) + cr(iage,ia,iab,is,ie,ir,1)  * fdis(iage,ia,iab,is,ie,ir,ii,1) /sum(muraceedu(:,:,ir))
                                do ip = 0, 1
                                    if (iage == 1) first_a(ir) = first_a(ir) +  a_bor(ia)  * fdis(1,ia,iab,is,ie,ir,ii,ip)
                                    cc       = cc       + cons(iage,ia,iab,is,ie,ir,ii,ip) * fdis(iage,ia,iab,is,ie,ir,ii,ip)  
                                    kknew    = kknew    + a_bor(ia)                        * fdis(iage,ia,iab,is,ie,ir,ii,ip)    
                                    totinctax= totinctax+ inctax(iage,ia,iab,is,ie,ir)     * fdis(iage,ia,iab,is,ie,ir,ii,ip)  
                                    tottrans = tottrans + trans(iage,ia,iab,is,ie,ir,ii,ip)* fdis(iage,ia,iab,is,ie,ir,ii,ip)  
                                    cratenew = cratenew + cr(iage,ia,iab,is,ie,ir,ip)      * fdis(iage,ia,iab,is,ie,ir,ii,ip)   
                                    iavg     = iavg     + allinc(iage,ia,iab,is,ie,ir)     * fdis(iage,ia,iab,is,ie,ir,ii,ip)   
                                    avgtaxinc= avgtaxinc+ alltaxinc(iage,ia,iab,is,ie,ir)  * fdis(iage,ia,iab,is,ie,ir,ii,ip)   
                                    ybar     = ybar     + y(iage,iab,is,ie,ir)             * fdis(iage,ia,iab,is,ie,ir,ii,ip)   
                                    num_cr(ir) = num_cr(ir) + cr(iage,ia,iab,is,ie,ir,ip)  * fdis(iage,ia,iab,is,ie,ir,ii,ip)  
                                    crrace(ir) = crrace(ir) + cr(iage,ia,iab,is,ie,ir,ip)  * fdis(iage,ia,iab,is,ie,ir,ii,ip)  /sum(muraceedu(:,:,ir) )
                                    aprime   = aprime   + (aplus(iage,ia,iab,is,ie,ir,ii,ip)) * fdis(iage,ia,iab,is,ie,ir,ii,ip)  * psi(iage+1,iab,ir)
                                    beq(ir)  = beq(ir)  + (aplus(iage,ia,iab,is,ie,ir,ii,ip)) * fdis(iage,ia,iab,is,ie,ir,ii,ip)  * (1-psi(iage+1,iab,ir))  
                                    if (iage <= nret) ylevel(iab,ir) = ylevel(iab,ir) + y(iage,iab,is,ie,ir) * fdis(iage,ia,iab,is,ie,ir,ii,ip)   / poplevel(iab,ir)
                                    if (iage >  nret) totalpen        = totalpen        + pen(iab,ir)       * fdis(iage,ia,iab,is,ie,ir,ii,ip)  
                                    if (a_bor(ia) <= 0.0d0) frac_debt(ir) = frac_debt(ir) + fdis(iage,ia,iab,is,ie,ir,ii,ip)   /sum(muraceedu(:,:,ir) )
                                    if (a_bor(ia) <= 0.0d0) tot_debt = tot_debt + fdis(iage,ia,iab,is,ie,ir,ii,ip) 
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    
    i = 1
    sort_inc = 0.0d0
    dis_relinc = 0.0d0
    do iage = 1, nage
        do iab = 1, niab
            do ie = 1, neta
                do ir = 1, nr
                    do ia = 0, na
                        do is = 1, ns
                            sort_inc(i) = alltaxinc(iage,ia,iab,is,ie,ir)
                            dis_relinc(i) = sum(fdis(iage,ia,iab,is,ie,ir,:,:)) 
                            i = i + 1
                        end do
                    end do
                end do
            end do
        end do
    end do

    call sort(sort_inc,order_relinc)
    dis_relinc = dis_relinc(order_relinc)

    cum_relinc = 0.0d0
    cum_relinc(1) = dis_relinc(1)

    do i = 2,size(dis_relinc)
        cum_relinc(i) = cum_relinc(i-1) + dis_relinc(i)
        if (cum_relinc(i) <= 0.2) then
            inc_cutoff(1) = sort_inc(i)
        else if (cum_relinc(i) > 0.2 .and. cum_relinc(i) <= 0.4) then
            inc_cutoff(2) = sort_inc(i) 
        else if (cum_relinc(i) > 0.4 .and. cum_relinc(i) <= 0.6) then
            inc_cutoff(3) = sort_inc(i) 
        else if (cum_relinc(i) > 0.6 .and. cum_relinc(i) <= 0.8) then
            inc_cutoff(4) = sort_inc(i)
        end if
    end do


    pi_v = 0.0d0
    scale_factor = 0.0d0
    do iage = 1, nage
        do ia = 0, na
            do iab = 1, niab
                do is = 1, ns
                    do ie = 1, neta
                        do ir = 1, nr
                            if (alltaxinc(iage,ia,iab,is,ie,ir) <= inc_cutoff(1)) then
                                pi_v(iage,ia,iab,is,ie,ir) = vic(1)
                            elseif (alltaxinc(iage,ia,iab,is,ie,ir) > inc_cutoff(1) .and. alltaxinc(iage,ia,iab,is,ie,ir) <= inc_cutoff(2)) then
                                pi_v(iage,ia,iab,is,ie,ir) = vic(2)
                            elseif (alltaxinc(iage,ia,iab,is,ie,ir) > inc_cutoff(2) .and. alltaxinc(iage,ia,iab,is,ie,ir) <= inc_cutoff(3)) then
                                pi_v(iage,ia,iab,is,ie,ir) = vic(3)
                            elseif (alltaxinc(iage,ia,iab,is,ie,ir) > inc_cutoff(3) .and. alltaxinc(iage,ia,iab,is,ie,ir) <= inc_cutoff(4)) then
                                pi_v(iage,ia,iab,is,ie,ir) = vic(4)
                            else
                                pi_v(iage,ia,iab,is,ie,ir) = vic(5)
                            end if
                            scale_factor = scale_factor +  pi_v(iage,ia,iab,is,ie,ir) * sum(fdis(iage,ia,iab,is,ie,ir,:,:))
                        end do
                    end do
                end do
            end do
        end do
    end do

    pi_v = pi_v * prorate/scale_factor


    
    do iab = 1, niab
        do ir = 1, nr
            pen(iab,ir) = penkappa*ylevel(iab,ir)
        end do
    end do



    stealinc = gamma*ybar  * prorate
    loseaux = 0.0d0
    do iage = 1, nage
        do ia = 0, na
            do iab = 1, niab
                do is = 1, ns
                    do ie = 1, neta
                        do ir = 1, nr
                            loseaux= loseaux + y(iage,iab,is,ie,ir)*pi_v(iage,ia,iab,is,ie,ir)*sum(fdis(iage,ia,iab,is,ie,ir,:,:))
                        end do
                    end do
                end do
            end do
        end do
    end do
    losefrac = stealinc/loseaux


    unben = 0.0d0
    llnew = 0.0d0
    taxbase = 0.0d0
    do iage = 1, nage
        do ia = 0, na
            do iab = 1, niab
                do ie = 1, neta
                    do ir = 1, nr
                        do ii = 1, 5
                            do ip = 0, 1
                                unben   = unben   + mpl * hbar *eff(iage,iab,ir) *eta(iab,ir,ie)* phi(ir) * fdis(iage,ia,iab,2,ie,ir,ii,ip)
                                llnew   = llnew   +       hbar *eff(iage,iab,ir) *eta(iab,ir,ie)          * fdis(iage,ia,iab,1,ie,ir,ii,ip)
                                taxbase = taxbase + mpl * hbar *eff(iage,iab,ir) *eta(iab,ir,ie)          * fdis(iage,ia,iab,1,ie,ir,ii,ip)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    expense = unben + cons_jail + tottrans
    revenue = sum(confis) + totinctax
    
    gg = revenue - expense
    taup = totalpen/taxbase
    
    out = kknew**(1.0-alfa)*(llnew**alfa)*bar ! total output
    
    kdev = abs(kk-kknew)/kk
    cdev = abs(prorate - proratenew)/prorate

    kk = (1-ww)*kk + ww*kknew
    ll = (1-ww)*ll + ww*llnew

    crate = (1-ww1)*crate + ww1*cratenew
    prorate = (1-ww1)*prorate + ww1*proratenew 

    inves =  (n_p+dep)*kknew
    exdem = (cc + inves + gg - out)/out

    write(*,*)'                               '
    write(*,'(A20, I4)')'Iteration times= ', iter
    write(*,'(2(A20, E10.3))')'Capital error= ', kdev, 'Crime rate error= ', cdev
    write(*,'(2(A20, E10.3))') 'Total error= ', max(kdev,cdev), 'Excess demand= ', exdem
    write(*,*)'__________________________________________________________________'
    write(*,fmt="(A20,3(A10))")' ','Moments', 'Target', 'Diff. (%)'
    write(*,fmt="(A20,3(F10.5))")'Property Crate%= ', prorate*100, 3.5461 , (prorate*100 - 3.5461)/ 3.5461 *100
    write(*,fmt="(A20,3(F10.5))")'Total Sto/GDP(%)= ', stealinc/out*100, 0.1398, (stealinc/out*100 - 0.1398)/0.1398 *100
    write(*,fmt="(A20,3(F10.5))")'% Black incar= ', incar_b/incar*100, 33.8, (incar_b/incar*100 - 33.8)/33.8 *100
    write(*,fmt="(A20,3(F10.5))")'% Black (a<=0)= ', frac_debt(1)*100, 19.25, (frac_debt(1)*100 - 19.25)/19.25 *100
    write(*,fmt="(A20,3(F10.5))")'% All (a<=0)= ', tot_debt*100, 10.36, (tot_debt*100 - 10.36)/10.36 *100
    write(*,*)'__________________________________________________________________'
    write(*,fmt="(A20,F10.5)")'Total Crate%= ',crate*100
    write(*,fmt="(A20,F10.5)")'Drug Crate%= ', drugrate*100
    write(*,fmt=*)'__________________________________________________________________'
    write(*,fmt="(3(A15, F7.2))")'20-29: ', sum(incargroup(1:2)), '40-54: ', sum(incargroup(5:7))
    write(*,fmt="(3(A15, F7.2))")'(Tar)20-29: ', 47.28, '40-54: ', 21.57
    write(*,fmt=*)'__________________________________________________________________'

    iter = iter + 1

    open(unit=10, file="../results/output_noscarring_drug.txt", status="old", action="write", position="append")
    write(unit=10,fmt=*)"         "
    write(unit=10,fmt=*)'Iter number= ',iter
    write(unit=10,fmt="(3(A20, F10.5))") 'K= ', kk, 'r= ', rbar, 'w= ', mpl
    write(unit=10,fmt="(3(A20, F10.5))") 'L= ', ll, 'C= ', cc, 'Output= ', out
    write(unit=10,fmt="(3(A20, F10.5))") '1+r= ', atr, 'Trans= ', tottrans, 'K/Y= ', kk/out
    write(unit=10,fmt="(3(A20, F10.5))") 'Trans cutoff= ', kappa*iavg, 'Avg total inc= ', iavg, 'Total Pen= ', totalpen
    write(unit=10,fmt="(4(A20, F10.5))") 'Inves= ', inves, 'Taup= ', taup, 'Tot inctax= ', totinctax
    write(unit=10,fmt="(4(A20, F10.5))") 'Beq(tot)= ', sum(beq),'Beq(b)= ', beq(1),'Beq(w)= ', beq(2)
    write(unit=10,fmt="(4(A20, F10.5))") 'Beqdis(b)= ', beqdis(1),'Beqdis(w)= ', beqdis(2) ,'Beqdis(b/w)= ', beqdis(1)/ beqdis(2)
    write(unit=10,fmt="(4(A20, F10.5))") '(b)= ', sum(muraceedu(:,:,1)),'(w)= ', sum(muraceedu(:,:,2)),'Avg Beq(b/w)= ', (beq(1)/sum(muraceedu(:,:,1)))/(beq(2)/sum(muraceedu(:,:,2)))
    write(unit=10,fmt="(4(A20, F10.5))") 'beqtrans(20,b)= ', beqtrans(1,1),'beqtrans(20,w)= ', beqtrans(1,2),'b/w= ', beqtrans(1,1)/beqtrans(1,2)
    write(unit=10,fmt="(4(A20, F10.5))") 'beqtrans(50,b)= ', beqtrans(31,1),'beqtrans(50,w)= ', beqtrans(31,2),'b/w= ', beqtrans(31,1)/beqtrans(31,2)
    write(unit=10,fmt="(3(A20, F10.5))") 'excess demand= ', exdem, 'Ybar= ', ybar, 'Crime rate= ', cratenew
    write(unit=10,fmt="(3(A20, F10.5))") 'kdev= ',kdev,'cdev= ',cdev, 'Error= ', max(kdev,cdev)
    close(unit=10)
    
end do


income = 0.0d0
preincome = 0.0d0
preincome(:,:,:,:,:,1) = y 
preincome(:,:,:,:,:,2) = y 
preincome(:,:,:,:,:,3) = y 
preincome(:,:,:,:,:,4) = 0.0d0
preincome(:,:,:,:,:,5) = y 
do iage = 1, nage
    do ia = 0, na
        do iab = 1, niab
            do is = 1,ns
                do ie = 1, neta
                    do ir = 1, nr
                        do ip = 0, 1
                            income(iage,ia,iab,is,ie,ir,1,ip) = y(iage,iab,is,ie,ir) + a_bor(ia)*rbar + trans(iage,ia,iab,is,ie,ir,1,ip)

                            income(iage,ia,iab,is,ie,ir,2,ip) = y(iage,iab,is,ie,ir) + a_bor(ia)*rbar + trans(iage,ia,iab,is,ie,ir,2,ip)

                            income(iage,ia,iab,is,ie,ir,3,ip) = y(iage,iab,is,ie,ir) + a_bor(ia)*rbar + trans(iage,ia,iab,is,ie,ir,3,ip)

                            income(iage,ia,iab,is,ie,ir,4,ip) = a_bor(ia)*rbar + trans(iage,ia,iab,is,ie,ir,4,ip)

                            income(iage,ia,iab,is,ie,ir,5,ip) = y(iage,iab,is,ie,ir) + a_bor(ia)*rbar + trans(iage,ia,iab,is,ie,ir,5,ip)
                        end do
                    end do
                end do
            end do
        end do
    end do
end do




open(unit=10, file="../results/output_noscarring_drug.txt", status="old", action="write", position="append")
write(unit=10,fmt=*)'         '
write(unit=10,fmt=*)'Final Results'
write(unit=10,fmt='(4(A20,F10.5))')'K= ',kknew,'Aprime= ',aprime,'eL= ',ll
write(unit=10,fmt='(4(A20,F10.5))')'Y= ',out,'C= ',cc,'mpl= ',mpl
write(unit=10,fmt='(4(A20,F10.5))') 'excess demand= ', exdem , 'losefrac= ', losefrac
write(unit=10,fmt='(4(A20,F10.5))') 'Crime rate= ',crate*100, 'Crime rate (B)= ', crrace(1)*100,'Crime rate (W)= ', crrace(2)*100
write(unit=10,fmt='(4(A20,F10.5))') 'Property Crate= ',prorate*100, 'Property Crate (B)= ', procrrace(1)*100,'Property Crate (W)= ', procrrace(2)*100
write(unit=10,fmt='(4(A20,F10.5))') 'Drug Crate= ',drugrate*100, 'Drug Crate (B)= ', drugcrrace(1)*100, 'Drug Crate (W)= ', drugcrrace(2)*100
write(unit=10,fmt='(4(A20,F10.5))')'Gov exp= ', gg,'Gov exp (% of Y)= ', gg/out*100, 'Total transfer= ', tottrans
write(unit=10,fmt='(4(A20,F10.5))')'Cost of cbar= ', cons_jail, 'Total incar.= ', incar
write(unit=10,fmt='(A50,3(F10.5))')'Total confiscated criminals labor earnings= ', confis(1) + confis(2)
write(unit=10,fmt='(A50,F10.5)')'Total confiscated stolen earnings= ', confis(3)

write(unit=10,fmt=*)'__________________________________________________________________'
write(unit=10,fmt="(A20,3(A10))")' ','Moments', 'Target ', 'Diff. (%)'
write(unit=10,fmt="(A20,3(F10.5))")'Property Crate%= ',prorate*100, 3.5461 , (prorate*100 - 3.5461)/ 3.5461 *100
write(unit=10,fmt="(A20,3(F10.5))")'Total Sto/GDP(%)= ', stealinc/out*100, 0.1398, (stealinc/out*100 - 0.1398)/0.1398 *100
write(unit=10,fmt="(A20,3(F10.5))")'% Black incar= ', incar_b/incar*100, 33.8, (incar_b/incar*100 - 33.8)/33.8 *100
write(unit=10,fmt="(A20,3(F10.5))")'% Black (a<=0)= ', frac_debt(1)*100, 19.25, (frac_debt(1)*100 - 19.25)/19.25 *100
write(unit=10,fmt="(A20,3(F10.5))")'% All (a<=0)= ', tot_debt*100, 10.36, (tot_debt*100 - 10.36)/10.36 *100
write(unit=10,fmt=*)'__________________________________________________________________'
write(unit=10,fmt="(A20,F10.5)")'Total Crate%= ',crate*100
write(unit=10,fmt="(A20,F10.5)")'Drug Crate%= ', drugrate*100
write(unit=10,fmt=*)'__________________________________________________________________'
write(unit=10,fmt="(3(A15, F7.2))")'20-29: ', sum(incargroup(1:2)), '40-54: ', sum(incargroup(5:7))
write(unit=10,fmt="(3(A15, F7.2))")'(Tar)20-29: ', 47.28, '40-54: ', 21.57
write(unit=10,fmt=*)'__________________________________________________________________'
close(unit=10)


open(unit=57, file = '../model_input/initial_noscarring_drug.txt', status="replace", action="write", position="rewind")
write(unit=57,fmt=*)kk    
write(unit=57,fmt=*)ll   
write(unit=57,fmt=*)crate 
write(unit=57,fmt=*)iavg
write(unit=57,fmt=*)taup 
write(unit=57,fmt=*)avgtaxinc
write(unit=57,fmt=*)losefrac
write(unit=57,fmt=*)ybar
write(unit=57,fmt=*)prorate
write(unit=57,fmt=*)beq(1)
write(unit=57,fmt=*)beq(2)
do iab = 1, niab
    do ir = 1, nr
        write(unit=57,fmt=*)pen(iab,ir)
    end do
end do
do iage = 1, nage
    do ia = 0, na
        do iab = 1, niab
            do is = 1, ns
                do ie = 1, neta
                    do ir = 1, nr
                        write(unit=57,fmt=*)pi_v(iage,ia,iab,is,ie,ir)
                    end do
                end do
            end do
        end do
    end do
end do
close(unit=57)
write(*,*)'Finish writing initial'



open(unit=33, file = '../results/simu_data_noscarring_drug', status="replace", action="write", position="rewind")

write(unit=33,fmt=*)nage
write(unit=33,fmt=*)na
write(unit=33,fmt=*)niab
write(unit=33,fmt=*)neta
write(unit=33,fmt=*)ns
write(unit=33,fmt=*)nr
write(unit=33,fmt=*)amax
write(unit=33,fmt=*)amin
write(unit=33,fmt=*)agrow
write(unit=33,fmt=*)rbar


do ia = 0, na
    write(unit=33,fmt=*) a_bor(ia)
end do



do iage = 1,nage
    do iab = 1, niab
        do is = 1, ns
            do ie = 1, neta
                do ir = 1, nr
                    do ia = 0, na
                        do ii = 1, 5
                            do ip = 0, 1
                                write(unit=33,fmt=*)fdis(iage,ia,iab,is,ie,ir,ii,ip)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do



close(unit=33)
write(*,*)'Finish writing simulation data'



call system_clock(endtime)
seconds = real(endtime - starttime)/real(rate)
write(*,fmt='(A30,I5)')"elapsed time(mins)= ", int(seconds/60.0d0)
write(*,fmt='(A30,I5)')"elapsed time(seconds)= ", int(mod(seconds,60.0d0))


contains


function saving(x)
    implicit none
    real*8, intent(in) :: x
    real*8 :: saving, ap, varphi
    

    ap = x

    if (ii == 1) then
        cons_com = max(dinc+atr*a_bor(ia)-ap+trans(iage,ia,iab,is,ie,ir,ii,ip)&
                    - inctax(iage,ia,iab,is,ie,ir)+ beqtrans(iage,ir), 1d-10)
    else if(ii == 2) then
        cons_com = max((1.0d0-losefrac)*dinc+atr*a_bor(ia)-ap +trans(iage,ia,iab,is,ie,ir,ii,ip)&
                    - inctax(iage,ia,iab,is,ie,ir)+ beqtrans(iage,ir), 1d-10)
    else if(ii == 3) then
        cons_com = max(dinc+atr*a_bor(ia)-ap +trans(iage,ia,iab,is,ie,ir,ii,ip)+steal(iage)&
                    - inctax(iage,ia,iab,is,ie,ir)+ beqtrans(iage,ir), 1d-10)
    else if(ii == 5) then
        cons_com = max((1.0d0-losefrac)*dinc+atr*a_bor(ia)-ap +trans(iage,ia,iab,is,ie,ir,ii,ip)+steal(iage)&
                        - inctax(iage,ia,iab,is,ie,ir)+ beqtrans(iage,ir), 1d-10)
    end if

    call linint_Grow(ap, amin, amax, agrow, na, ial, iar, varphi)
    ial = max(min(ial, na-1), 0)
    iar = max(min(iar, na), 1)
    varphi = max(min(varphi, 1.0d0), 0.0d0)

    saving = 0.0d0

    saving = (varphi*ev(iage+1,ial,iab,is,ie,ir) + (1-varphi)*ev(iage+1,iar,iab,is,ie,ir))**(1-sigma1)/(1-sigma1)
    
    saving = - (cons_com**(1.0d0-sigma1)/(1.0d0-sigma1) + beta * psi(iage+1,iab,ir) * saving)

end function



subroutine interpolate(ij)
    implicit none
    integer, intent(in) :: ij

    do ia = 0, na
        do iab = 1, niab
            do is = 1, ns
                do ie = 1, neta
                    do ir = 1, nr
                        ev(ij,ia,iab,is,ie,ir) = 0.0d0
                        do is_p = 1, ns
                            do ie_p = 1, neta
                                
                                ev(ij,ia,iab,is,ie,ir) = ev(ij,ia,iab,is,ie,ir) + pis(ij,is,is_p,ir)*pieta(iab,ir,ie,ie_p)*&
                                (drugpro(ij,ir)*val(ij,ia,iab,is_p,ie_p,ir,1) + (1-drugpro(ij,ir))*val(ij,ia,iab,is_p,ie_p,ir,0))
                            end do
                        end do
                        
                        ev(ij,ia,iab,is,ie,ir) = (ev(ij,ia,iab,is,ie,ir)*(1-sigma1))**(1/(1-sigma1))
                    end do
                end do
            end do
        end do
    end do

end subroutine

subroutine transfer(x, incwel)
    implicit none
    real*8, intent(in) :: incwel
    real*8, intent(out) :: x

    if (incwel > kappa * iavg) then
        x = 0.0d0
    else
        x  = kappa * iavg - incwel
    end if

end subroutine


subroutine incometax(tottaxinc, paidtax)
    implicit none
    real*8, intent(in) :: tottaxinc
    real*8, intent(out) :: paidtax
    real*8 :: relinc

    if (tottaxinc > 1d-10) then
        relinc = tottaxinc/avgtaxinc
        paidtax = (1-tau1*relinc**(-tau2))*tottaxinc
    else
        paidtax = 0.0d0
    end if

end subroutine



end program racialgap























