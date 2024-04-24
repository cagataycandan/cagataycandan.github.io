%% EE 306 : An Application of Bayes Theorem in Communications
% We revisit the problem discussed in the lecture hours. 
%
% Assume that _S_ is a binary random variable taking values 
% 0 and 1 with the probability of _p_ and _q_ = 1 - _p_, respectively. 
%
% $$ S = \left\{ \begin{array}{lll} 1,  & \ & \textrm{with probability p} \\ 0, & \ &\textrm{with probability 1-p = q} \end{array} \right. $$
%
% We observe not _S_, but its noise corrupted version:
% 
% $$ X_k = S + N_k \quad k=\{1,2, \ldots, N \} $$ 
% 
% Here $X_k$ is the k'th observation. $N_k$ is the random variable
% representing noise which corrupts the signal _S_. We assume that the
% random variable $N_k$ is independent and identically Gaussian distributed
% with zero mean and variance $\sigma_n^2$,
% 
% $$ N_k \sim N_{n_k}(0, \sigma_n^2)  = \frac{1}{\sqrt{2\pi}} 
% \exp \left( -\frac{n_k^2}{2 \sigma_n^2} \right). $$
%
% One can visualize the situation as follows. The transmitter sends either
% 1 or 0; which is the information bit to be delivered to the receiver. The
% receiver observes the signal _S_ in the presence of noise. (The noise is
% Gaussian distributed as noted above.) Assume that $X_1$ is 1.01. If $S$
% is 1, the noise $N_1$ has to be 0.01 to produce this observation.
% Similarly, if $S$ is 0, the noise $N_1$ has to be 1.01 to produce the
% same observation. Our goal is to decide on S given an observation, say
% $X_1$. It is intuitively clear that for $X_1 = 1.01$, the event of $S=1$
% is more likely than $S=0$. We would like to quantify this intuition by
% working out the probability of $S$ given $X_1$, in terms of density
% functions it is $f_{S|X_1}(s|x_1)$. 
% 
% The distribution of the $S$ before we observe $X_1$ is called the
% a-priori distribution of $S$. From the problem statement, the a-priori
% density can be written as
% 
% $$ f_S(s) = p_0\delta(s) + p_1\delta(s-1). $$
%
% The conditional density $f_{S|X_1}(s|x_1)$ can be written as 
% 
% $$ f_{S|X_1 }(s|x_1) = \frac { f_{X_1|S}(x_1 | s) f_S(s) } 
%    { f_{X_1} (x_1) }. $$ 
%
% The density on the left, $f_{S|X_1 }(s|x_1)$ is called the posterior density. 
% Stated differently, our goal is to calculate the posterior density. 
% 
% From the observation model, the random variable $X_1$ given $S=s$ is distributed as
% 
% $$ X_1 | s  \sim N_{x_1}(s, \sigma_n^2) $$ 
%
% This result is easy to see, once we interpret this result as fixing the random variable _S_ to a
% constant numerical value shown as _s_. In other words, $X_1$ given $S=s$ is a new
% random variable and _s_ is just a parameter of this distribution (it
% is not any more a random variable after fixing its value).  
%
% Then, the posterior density can be written as $f_{S|X_1 }(s|x_1)$
%
% $$f_{S|X_1 }(s|x_1) = \frac{ N_{x_1}(s, \sigma_n^2) (p_0\delta(s) +
% p_1\delta(s-1) ) }  
%    { f_{X_1} (x_1) } 
% = 
% \frac{ p_0 N_{x_1}(0, \sigma_n^2) \delta(s)  + p_1 N_{x_1}(1, \sigma_n^2)\delta(s-1) }
%    { f_{X_1} (x_1) } 
% $$
%
% The denominator of the ratio given above can be calculated from the relation
%
% $$ f_{X_1} (x_1) = \int_{-\infty}^{\infty} f_{X_1 , S }(x_1, s) ds = 
% \int_{-\infty}^{\infty} f_{X_1 | S }(x_1|s) f_S(s)
% ds $$ 
% 
% which is simply the integral of the numerator of the posterior
% density. 
% 
% $$ f_{X_1} (x_1) = \int_{-\infty}^{\infty} 
% \left(p_0 N_{x_1}(0, \sigma_n^2) \delta(s)  + p_1 N_{x_1}(1, \sigma_n^2)\delta(s-1)
% \right) ds
% = p_0 N_{x_1}(0, \sigma_n^2) + p_1 N_{x_1}(1, \sigma_n^2) 
% $$ 
%
% It should be noted that the denominator is just a scalar (a constant
% value) which does not functionally affect the density on the right hand
% side; but it is present to scale the right hand side of the equation such
% that the the area under the density is normalized to 1. It can be immediately 
% verified that with this normalization the area under the density is indeed 1. 
%
% The posterior density can be finalized as 
% 
% $$ f_{S|X_1 }(s|x_1) = \frac{ p_0 N_{x_1}(0, \sigma_n^2)}   
%  { p_0 N_{x_1}(0, \sigma_n^2) + p_1 N_{x_1}(1, \sigma_n^2) } \delta(s) 
% + 
% \frac{ p_1 N_{x_1}(1, \sigma_n^2)}   
%  { p_0 N_{x_1}(0, \sigma_n^2) + p_1 N_{x_1}(1, \sigma_n^2) } \delta(s-1). 
% $$
% 
% Note that, the posterior distribution updates the probability of the event $S=0$ from
% $p_0$ to $\widehat{p}_0$ by processing the observation. Similarly, $p_1$ is updated 
% to $\widehat{p}_1$:
% 
% $$ \widehat{p}_0 = \frac{ p_0 N_{x_1}(0, \sigma_n^2)}   
%  { p_0 N_{x_1}(0, \sigma_n^2) + p_1 N_{x_1}(1, \sigma_n^2) }, 
% $$
% 
% $$ \widehat{p}_1 = \frac{ p_1 N_{x_1}(1, \sigma_n^2)}   
%  { p_0 N_{x_1}(0, \sigma_n^2) + p_1 N_{x_1}(1, \sigma_n^2) }, 
% $$
%
% Hence, the prior density is updated from 
%
% $$ f_S(s) = p_0\delta(s) + p_1\delta(s-1), $$ 
%
% to the posterior
% 
% $$ f_{S|X_1}(s|x_1) = \widehat{p}_0 \delta(s) + \widehat{p}_1
% \delta(s-1). $$ 
% 
%
% Now, we numerically study this mentioned update using MATLAB's
% computational facilities. 
% 
% To do that, we first define the Gaussian distribution
% 
pdfNormal = @(x,mu,var) 1/sqrt(2*pi*var)*exp(-(x-mu).^2/2/var); 

%% 
% Let's do a sanity check and calculate the area under the pdf. To do that
% we fix the mean value to 0, variance to 1 and calculate the area under pdf with the quad
% function as follows: 
% 
area = quad(@(x)pdfNormal(x,0,1),-5,5), 

%%
% We only calculate the area under 5 standard deviations around the mean. It seems that everything is in order!
% 
% Let's define, the update equations:
update_p0 = @(p0,x1,sigma_n_sq) p0*pdfNormal(x1,0,sigma_n_sq)/... 
    (p0*pdfNormal(x1,0,sigma_n_sq) + (1-p0)*pdfNormal(x1,1,sigma_n_sq)); 
update_p1 = @(p1,x1,sigma_n_sq) 1 - update_p0(1-p1,x1,sigma_n_sq); 
%
%%
%
% Let's take $\sigma_n^2 = 2$, and p0 = p1 = 1/2 (Hence, the events of S=0 or S=1
% are equally likely, which is the a-priori information.)
sigma_n_sq = 2; 
p0 = 1/2; 
p1 = 1/2; 
%%
% 
%
% For now, assume that the observation $x_1$ is 1.01. For this observation,
% the updated probabilities can be calculated as 
% 
x1 = 1.01; 
[update_p0(p0,x1,sigma_n_sq) update_p1(p1,x1,sigma_n_sq)]

%%
% Given, the observation of x1=1.01; the updated probability for 
% s=1 is increased from its a-priori value of 0.5 to 0.5634. The increase is not significant for this
% case. This is due to large noise variance. That is, a there exists a
% fairly large probability for noise to be 1.01; therefore the event of S=0
% has also a fairly large probability. 
% 
% Let's decrease the noise variance and repeat the same experiment: 
sigma_n_sq = 1;
x1 = 1.01; 
[update_p0(p0,x1,sigma_n_sq) update_p1(p1,x1,sigma_n_sq)]

%%
sigma_n_sq = 1/2;
x1 = 1.01; 
[update_p0(p0,x1,sigma_n_sq) update_p1(p1,x1,sigma_n_sq)]

%%
sigma_n_sq = 1/4;
x1 = 1.01; 
[update_p0(p0,x1,sigma_n_sq) update_p1(p1,x1,sigma_n_sq)]

%%
% Note that, with 95% percent probability, a random pick from the Gaussian distribution
% lies in two standard deviations of the mean value. For the last
% case of $\sigma_n^2 = 1/4$, the two standard deviation interval around the
% mean is [-1,1]. Hence, 95% percent of time, we can not see the
% observation of $x_1 = 1.01$ when $s=0$. 
%
% In short, as the noise variance decreases, the probability of
% having a large valued noise also decreases and 
% observing $x_1=1.01$ for the event of $s=0$ becomes less and less likely. 
%
%
%
% Let's use a high noise variance and make multiple observations on the
% unknown symbol S. Let's take $\sigma_n^2 = 4$ and generate 5 noisy
% observations on $S$. 
%
% Without any loss of generality, let's examine the case of $s=0$

sigma_n_sq = 4;
s = 0; 
x = zeros(1,5); 
for k=1:5,
    x(k)  = s + sqrt(sigma_n_sq)*randn(1);
end;
x, 
%%
% It can be noted that the observations are large valued due to large noise
% variance. 
% 
% We would like to estimate the _s_ value given all observations. In other
% words, we would like to evaluate 
% 
% $$ f_{S|X_1,X_2, ..., X_5}(s|x_1, x_2, ... , x_5) $$
% 
% To evaluate the posterior density with 5 observations, we can do the
% updates in a recursive fashion. That is, we first find 
%
% $$ f_{S|X_1}(s|x_1) $$
% 
% then treat the updated probabilities of $\widehat{p}_0$ and
% $\widehat{p}_1$ as the prior probabilities before the observation of
% $X_2$
% and generate $f_{S|X_1,X_2}(s|x_1,x_2)$. 
%
% $$ f_{S|X_1,X_2}(s|x_1,x_2) = \frac{1}{A}f_{X_2|S,X_1}(x_2 | s, x_1) f_{S|X_1}(s|x_1) $$
%
% where $A$ is a constant normalizing the density to unit area, as shown above.  
% Note that the knowledge of $s$ and $x_1$ is equivalent to the knowledge
% of $s$ and $n_1$. The knowledge of $n_1$ (noise of the first observation)
% does not help at any other observations, since noise is independent at
% every observation. Hence, we can discard the conditioning event of $x_1$
% in this relation, since it does not affect the density on the right
% hand side:
%
% $$ f_{X_2|S,X_1}(x_2 | s, x_1) = f_{X_2|S}(x_2 | s ). $$
%
% This fact is stated as $X_2$ is _conditionally indendepent_ of $X_1$ given
% $S$. 
%
% With this observation, we can write the posterior density as 
% 
% $$ f_{S|X_1,X_2}(s|x_1,x_2) = \frac{1}{A}f_{X_2|S}(x_2 | s) f_{S|X_1}(s|x_1) $$
%
% In the last equation $f_{S|X_1}(s|x_1)$ is the updated density after the
% observation of $x_1$ which is the a-priori information before processing
% the second observation $x_2$. The other term $f_{X_2|S}(x_2 | s)$ is
% identical to the term appearing in the update equation for the single
% observation case. (This term is called the _likelihood_ term)
% 
% Hence, the only change between single observation and multiple
% observations case is taking the aposteriori distribution of the earlier
% iteration as the apriori density of the next iteration. 
%
% Given this discussion, we can calculate the posterior probability after
% each observation as: 
%
p0old = 1/2; p1old = 1/2; 
p0_posterior = zeros(1,5); p1_posterior = zeros(1,5); 
for k=1:5, 
    p0_posterior(k) = update_p0(p0old,x(k),sigma_n_sq); 
    p1_posterior(k) = update_p1(p1old,x(k),sigma_n_sq);
    p0old = p0_posterior(k);
    p1old = p1_posterior(k);
end; 
p0_posterior,

%%
% We can note that as more and more observations are collected, the
% probability of making correct decision shows a tendency to increase. 
%
%
% Let's repeat the same example with 30 observations
sigma_n_sq = 4;
s = 0; 
x = zeros(1,30); 
for k=1:30,
    x(k)  = s + sqrt(sigma_n_sq)*randn(1);
end;
x,
%%
% Let's do the posterior calculation:

p0old = 1/2; p1old = 1/2; 
p0_posterior = zeros(1,30); p1_posterior = zeros(1,30); 
for k=1:30, 
    p0_posterior(k) = update_p0(p0old,x(k),sigma_n_sq); 
    p1_posterior(k) = update_p1(p1old,x(k),sigma_n_sq);
    p0old = p0_posterior(k);
    p1old = p1_posterior(k);
end; 
p0_posterior,

%%
% Please compare the entries of the observation vector $x$ and the
% posterior density vector and try to understand why for some observations, there is 
% a significant increase in the probability immediately after the update. 
%
% The main goal of communication system design is to establish a reliable
% communication between parties with a minimal repetition of the
% transmitted symbols, that is achieving a reliable communication (low probability of error) 
% at a highest rate possible (with a low redundancy). 
% 
