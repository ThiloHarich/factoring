package factoring.fermat.lehman;

import java.util.Collection;

import factoring.fermat.FermatFact;

/****
 Code by Warren D. Smith, Dec 2011, to implement Lehman integer-factorizationByFactors
 algorithm.  Rigorous O(N^(1/3)) step factorizationByFactors (also prime-proving) algorithm
 invented by RS Lehman 1974.
 COMPILE:
 ========
 gcc LehmanClean.c -Wall -O6 -o lehcl
 gcc LehmanClean.c -Wall -DNDEBUG -O6 -o lehcl
 SPEED:
 ======
 79 microsec avg factor time for "hard" 43-bit integers (zero failures in
 1 million factorizations tested) using LehTune=3.0 and HartOLF=FALSE (DoTrial=FALSE)
 on Apple iMac 2.0 GHz intel core duo.  If DoTrial=TRUE then 190 microsec avg factor time.
 In contrast, a (poor) implementation of SQUFOF that I stole from Dario Alpern plus
 then improved to speed up its inner loop, performs the same 1 million factorizations
 in 82 microsec avg factor time but with a 4.5% failure rate.  I think the best SQUFOF
 implementations run about 1.5X faster and fail way less often than mine.  One way
 to reduce SQUFOF's failure rate is to add a trial division (TD) stage to it; but mine
 has no trial division.
 "Hard" Ns are products of two random primes with half the number of bits each.
 Average factorizationByFactors times to factor "hard" Ns:
 #bits..LEHMAN(TD).......LEHMAN(No TD)....PoorSQUFOF
 43.....190(tune=2.8)....85.9(tune=3.0)...82usec(4.5% fail)
 42.....152..............56.8(tune=4.0)...69....(4.6% fail)
 41.....122..............41.5(tune=5.0)...58....(4.6% fail)
 40......98..............31.8(tune=6.0)...49....(4.7% fail)
 39......79..............24.3(tune=6.0)...41.5..(4.8% fail)
 38......63..............18.7(tune=6.0)...35....(4.9% fail)
 37......51(tune=2.8)....14.2(tune=6.0)...30....(4.9% fail)
 All Lehman runs shown have zero failures in 1 million factorizations.
 (A "failure" is failing to return a valid nontrivial factor of a composite N,
 or failing to return N when N is prime, though all our N in this test are composite.
 Lehman *without* TD *would* sometimes fail if given "easy" N, but all N in this test hard.)
 Above Lehman runs with ASSERTs turned on; turning them off will speed up Lehman 3-11%, e.g.
 the following are (re)done with asserts off:
 #bits..LEHMAN(full TD)..LEHMAN(No TD)....PoorSQUFOF
 43.....182(tune=2.8)....78.9(tune=3.0)...82usec(4.5% fail)
 37......49.4(tune=2.8)..13.9(tune=6.0)...30....(4.9% fail)
 36......39.8(tune=2.8)..10.5(tune=6.0)...25.4..(5.0% fail)
 35......32.1(tune=2.8)..7.26(tune=7.0)...21.6..(5.1% fail)
 34......25.7(tune=2.65)..................18.6..(5.3% fail)
 33......20.8(tune=2.8)...................16.1..(5.3% fail)
 32......16.8(tune=2.8)...................13.9..(5.6% fail)
 31......13.48(tune=2.65).................11.7..(5.6% fail)
 30......10.91(tune=2.65)..................9.22.(5.6% fail)
 29......8.88(tune=2.65)...................7.32.(6.2% fail)
 28......7.15(tune=2.65)...................6.13.(6.7% fail)
 EFFECT OF W.B.HART'S (APPARENTLY USELESS) OLF SPEEDUP HEURISTIC:
 =================================================================
 All above runs with with HartOLF=0.  If we turn it on, here is what happens:
 32-bit hard N, asserts off, Tune=2.75, TD=TRUE:
 HartOLF:  0      0.001   0.01   0.1    0.2    0.4     0.8
 AvgTime:  16.73  16.72  17.23  18.69  19.15  19.54   19.87 usec.
 40-bit hard N, asserts off, Tune=2.75, TD=FALSE:
 HartOLF:  0    0.01    0.1    0.2   0.4   0.8
 AvgTime: 41.52 42.17  51.31  53.64 55.13 57.04
 As you can see, directly contradicting Hart's paper's claim that HartOLF will
 speed up Lehman dramatically (he says as much as 500%),
 I actually find that HartOLF always makes it slower, and the
 more OLFing you do, the worse it gets, monotonically.  This is for "hard N."  I've
 also tested other N-distributions; and with the others, HartOLF looks even worse than
 it does here.  (With old versions of this code I was occasionally able to find some
 locales where HartOLF helped slightly, but with the current version I know of no such locale.)
 HUGE WIN USING CutFrac=0.1:
 ===========================
 All above runs were using CutFrac=1.0 and are disappointing if DoTrial=TRUE.
 However, it then occurred to me that we could employ CutFrac=0.1 to get a
 speedup.  Using 0.1, Lehman will omit the last 90% of the trial factorizationByFactors then just go directly
 to square-finding.  If square finding fails (which is rare), Lehman then returns to
 trial dividing to do the missing 90%.   Here is what happens:
 40-bit hard composite N, asserts off, Tune=2.75:
 CutFrac:  1.0     0.1   0.05   NoTrialDivision   RECOMMENDED(0.1)
 Time:     93.7   48.27  45.12         41.4             48.27
 40-bit random N not divisible by 2 or 3 (include many prime N), asserts off, Tune=2.75:
 CutFrac:   1.0   0.4    0.2    0.1    0.01   RECOMMENDED(0.1)
 Time:     24.66 25.16  26.29  27.91   35.06        27.91
 40-bit N=A*B, where A,B 20-bit randoms not divisible by 2 or 3, asserts off, Tune=2.75:
 FirstCut:  1.0    0.8    0.6    0.4    0.2   0.1    0.01   RECOMMENDED(0.1)
 Time:     28.22  27.99  27.81  27.85  28.61  29.95  38.26       29.95
 For most serious uses of Lehman you are going to apply a pseudoprimality test before
 trying to factor N so that N will almost always be composite.  I therefore
 recommend CutFrac=0.1 for general purpose use.
 With CutFrac=0.1, Lehman is clearly superior to my SQUFOF for 42-bit integers (or less).
 (For 43 bits, it may be debatable which is superior.)
 For Ben Buhrman's best-in-world SQUFOF implementation, the crossover point also
 appears to be 43 bits provided optimum tuning and CutFrac both are employed by Lehman.
 (However that claim is preliminary and Buhrman may do more testing later.)
 TESTING:
 ========
 (Some of these comments about tests may be out of date since pertain to older code versions,
 but hopefully still true.)
 Experimentally: Always returns N if N=prime<65536 and DoTrial=FALSE and 0.1<=Tune<=9.6.
 But: If Tune<=0.8 then the   gcd(a+b, N)   can return N as its mode of success if N=3.
 I.e. the theorem that the special prime-detect return always will fire for N prime can thus
 be violated if Tune<0.8 and N=3 (but still correctly returns N in these cases), but it's
 valid if Tune>=0.8  and N>2.
 Also valid if Tune>=0.25 and N>31.
 Also valid if Tune>=0.10 and N>83.
 Current code forces trial division for primes<=83 unless explicitly turned off.
 Using Tune too large will cause error exits galore due to overflows.  Do not use
 Tune>9, that regime seems useless at best and bogus at worst.   Tune=2.8 to 9.0
 can cause overflow if N is large enough.
 High end: if N < 2^43.000049 = 8796393022207  then with good choice of Tune
 integer overflows are extremely rare.
 But I nevertheless have seen such rare overflows with some code versions and some parameters
 for N as small as  1991049599221 = 2^40.85
 but I currently believe that the rare overflows that arise for N<2^43.01 do not
 ever hurt correctness.  (Comments inside code discuss why I believe that.)
 Factor Validity: testing 100 million N with N=A*B with
 2611953 <= A <= B <= 2965920
 LehmanFactor always succeeded in factorizationByFactors.
 This is highly reassuring but could be deceptive since
 (a) changing the tuning constants could alter the failure set.
 (b) It appears that Lehman can be sensitive to the distribution of N and
 there are "hard locales" and "easy locales" which are not obvious.
 ***/
public class LehmanYafuFact extends FermatFact {

	static boolean[] issq1024;
	static byte[] issq4199;
	private double tune;

	public LehmanYafuFact(float tune) {
		this.tune = tune;
	}

	static long gcd64(long a, long b){
		long t;
		if(a<b){ t=b; b=a; a=t; }
		while(b != 0){
			t = b;
			b = a % b;
			a = t;
		}
		return( a );
	}

	static void MakeIssq(){
		issq1024 = new boolean[1024];
		issq4199 = new byte [4199];
		int i;
		for(i=0; i<1024; i++){ issq1024[i] = false; }
		for(i=0; i<1024; i++){ issq1024[(i*i)%1024] = true; }
		for(i=0; i<4199; i++){ issq4199[i] = 0; }
		for(i=0; i<3465; i++){ issq4199[(i*i)%3465] |= 2; }
		for(i=0; i<4199; i++){ issq4199[(i*i)%4199] |= 1; }
		System.out.println("Square tables built.\n");
	}

	static int[] prime = new int [6543]; //the 6542 primes up to 65536=2^16, then sentinel 65535 at end

	static void MakePrimeTable(){
		int i,j,k;
		prime[0]=2;
		prime[1]=3;
		prime[2]=5;
		k=3;
		for(i=7; i<65536; i+=2){
			boolean isPime = true;
			for(j=0; prime[j]* prime[j] <= i && isPime; j++){
				if(i%prime[j]==0)
					isPime = false;
			}
			if (isPime) {
				prime[k] = i;
				k++;
			}
		}
		assert(k==6542);
		prime[k] = 65535; //sentinel
		System.out.printf("Prime table[0..%d] built: ", k);
		for(i=0; i<20; i++){ System.out.printf("%d,", prime[i]); }
		System.out.printf("%d,...,%d,(%d)\n", prime[20],prime[6541],prime[6542]);
	}

	@Override
	public long findFactors(long n, Collection<Long> primeFactors) {
		if (issq1024 == null) {
			MakeIssq();
			MakePrimeTable();
		}
		final long factor = LehmanFactor(n, tune, 0.0, true, 1, primeFactors);
		//        long factor = LehmanFactor(n, 2.5, 0.0, true, 1.0);
		if (factor != n && primeFactors != null) {
			primeFactors.add(n / factor);
		}
		return factor;
	}

	static long LehmanFactor(long N, double Tune, double HartOLF, boolean DoTrial, double CutFrac, Collection<Long> primeFactors){
		int b,p,k,r,B,U,Bred,Bred2,inc,FirstCut,ip;
		long a,c,kN,kN4,B2,N480,UU;
		double Tune2, Tune3, x;
		if((N&1)==0) return(2); //N is even
		if(Tune<0.1){
			System.out.printf("Sorry0, Lehman only implemented for Tune>=0.1\n");
			return(0);
		}
		B = (int) (Tune * Math.pow(N, 1.0/3.0));
		FirstCut = (int) (CutFrac*B);
		if(FirstCut<84){ FirstCut=84; } //assures prime N will not activate "wrong" Lehman return
		if(FirstCut>65535){ FirstCut = 65535; }

		if(DoTrial){
			for(ip=1; ; ip++){ //trial division
				p = prime[ip];
				if(p>FirstCut) break;
				if(N%p==0) {
					if (primeFactors == null || N == p)
						return p;
					primeFactors.add((long) p);
					N = N/p;
				}
			}
		}
		B = (int) (Tune * Math.pow(N, 1.0/3.0));

		if(N>=8796393022207L){
			System.out.printf("Sorry1, Lehman only implemented for N<8796393022207\n");
			return(0);
		}
		Tune2 = Tune*Tune;
		Tune3 = Tune2*Tune;
		Bred = (int) (B / Tune3);

		if(HartOLF>0.0){ // Hart's "OLF" algorithm is tried to get more speed on average...?
			assert(gcd64(480, N)==1);
			N480 = N*480;
			UU = Long.MAX_VALUE / N480;
			Bred2 = (int) (B * HartOLF);
			if(UU>Bred2) UU=Bred2;
			UU *= N480;
			for(kN=N480; kN<UU; kN+=N480){
				a = (long) (Math.sqrt(kN)+0.999999999999999);
				c = a*a-kN;
				assert(c>=0);
				if(issq1024[(int)(c & 1023)]){
					if((issq4199[(int)(c%3465)]&2) != 0){
						if((issq4199[(int)(c%4199)]&1) != 0){
							b = (int) Math.sqrt(c+0.9);
							if(b*b==c){
								if(a>=b){
									B2 = gcd64(a-b, N);
									if(1<B2 && B2<N){ return(B2); }
								}
								B2 = gcd64(a+b, N);
								if(1<B2 && B2<N){ return(B2); }
							}
						}
					}
				}
				//				if(PrimeMath.isSquare(c)){
				//					b = (int) Math.sqrt(c+0.9);
				//					if(b*b==c){
				//						if(a>=b){
				//							B2 = gcd64(a-b, N);
				//							if(1<B2 && B2<N){ return(B2); }
				//						}
				//						B2 = gcd64(a+b, N);
				//						if(1<B2 && B2<N){ return(B2); }
				//					}
				//				}
			}
		}

		B2 = B*B;
		kN = 0;

		//Lehman suggested (to get more average speed) trying highly-divisible k first. However,
		//my experiments on trying to to that have usually slowed things down versus this simple loop:
		for(k=1; k<=Bred; k++){
			if((k&1) != 0)
			{ inc=4; r= (int) ((k+N)%4); }
			else{ inc=2; r=1; }
			kN += N;
			assert(kN == k*N);
			if(kN >= Long.MAX_VALUE){
				System.out.printf("Sorry2, overflow, N=%llu is too large\n", N);
				return(0);
			}
			//Actually , even if overflow occurs here, one could still use approximate
			//arithmetic to compute kN4's most-signif 64 bits only, then still exactly compute x...
			//With appropriate code alterations is should be possible to extent the range... but
			//I have not tried that idea for trying to "cheat" to gain more precision.
			kN4 = kN*4;
			x = Math.sqrt(kN);
			a = (long) x;
			if(a*a==kN){
				B2 = gcd64(a, N);
				assert(B2>1);
				assert(B2<N);
				return(B2);
			}
			x *= 2;
			a = (long) (x+0.9999999665); //very carefully chosen.
			//            a = (long) Math.ceil(x); //very carefully chosen.
			//Let me repeat that: a = x+0.9999999665.  Really.
			b= (int) (a%inc);  b = (int) (a + (inc+r-b)%inc);   //b is a but adjusted upward to make b%inc=r.
			assert( b%inc == r );
			assert( b>=a );
			assert( b<=a+4 );
			c = (long)b*(long)b - kN4;  //this is the precision bottleneck.
			//At this point, I used to do a test:
			//if( c+kN4 != (long)b*(long)b ) //overflow-caused failure: exit!
			//	System.out.printf("Sorry3, unrepairable overflow, N=%llu is too large\n", N);
			//  return(0);
			//However, I've now reconsidered.  I claim C language computes c mod 2^64 correctly.
			//If overflow happens, this is irrelevant because c is way smaller than 2^64, kN4, and b*b.
			//Hence c should be correct, despite the overflow. Hence I am removing this error-exit.

			U = (int) (x + B2/(2*x));
			//old code was  U = SQRT((double)(B2+kN4+0.99));   and was 4% slower.

			//Below loop is: for(all integers a with 0<=a*a-kN4<=B*B and with a%inc==r)
			for(a=b;  a<=U;  c+=inc*(a+a+inc), a+=inc ){
				//again, even though this assert can fail due to overflow, that overflow should not matter:
				//assert( c == (long)a*(long)a-kN4 );
				/** Programming trick:    naive code:     c = a*a-kN4;
                 In the inner loop c is bounded between 0 and T^2*N^(2/3)
                 and can be updated additively by adding inc*(anew+aold) to it
                 when we update a to anew=a+inc. This saves a multiplication in
                 the inner loop and/or allows us to reduce precision. **/
				if(issq1024[(int)(c & 1023)]){
					if((issq4199[(int)(c%3465)]&2) != 0){
						if((issq4199[(int)(c%4199)]&1) != 0){
							b = (int) Math.sqrt(c + 0.9);
							if(b*b==c){ //square found
								B2 = gcd64(a+b, N);
								assert(B2>1);
								//                                if(B2>=N){ System.out.printf("theorem failure: B2=%llu N=%llu\n", B2,N); }
								if (primeFactors == null)
									return(B2);
								if (Tune == 1) {
									primeFactors.add(B2);
									if (N != B2)
										primeFactors.add(N/B2);
									return 1;
								}
							}
						}
					}
				}
				//				if(PrimeMath.isSquare(c)){
				//					b = (int) Math.sqrt(c + 0.9);
				//					if(b*b==c) { //square found
				//						B2 = gcd64(a + b, N);
				//						assert (B2 > 1);
				//						if (B2 >= N) {
				//							System.out.printf("theorem failure: B2=%llu N=%llu\n", B2, N);
				//						}
				//						return (B2);
				//					}
				//			}
			}
		}

		//square-finding has failed so resume missing part of trial division:
		if(DoTrial){
			if(B>65535) B = 65535;
			for(ip=0 ; ; ip++){
				p = prime[ip];
				if(p>=B) break;
				if(N%p==0)
					return(p);
			}
		}

		return(N); //N is prime
	}


	public static void main(String [] args){
		long N, M;
		MakeIssq();
		MakePrimeTable();


		//Here are some typical calls to LehmanFactor.
		//  LehmanFactor(N, (tune from 0.1 to 9.6), (tune from 0 to 5.0),
		//            (TRUE unless want to skip trial factorizationByFactors which would be unusual),
		//            (TRUE if want to try OLF speculative speedup FALSE if skip it) );
		N=3141592651L;
		M = LehmanFactor(N, 2.5, 0.0, true, 0.4, null);
		System.out.println("A factor of " + N + " is " + M);

		N=3141592661L; //prime
		M = LehmanFactor(N, 2.5, 0.0, true, 0.5, null);
		System.out.println("A factor of " + N + " is " + M);

		N = 7919; N *= 10861;
		M = LehmanFactor(N, 1.0, 0.0, true, 0.1, null);
		System.out.println("A factor of " + N + " is " + M);

		N =  1299709; N *=  2750159;
		M = LehmanFactor(N, 1.0, 0.0, false, 0.1, null);
		System.out.println("A factor of " + N + " is " + M);
	}
}
