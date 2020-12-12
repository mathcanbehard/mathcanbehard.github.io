---
layout: post
title:  "Two confidence intervals for the rate parameter of the Exponential distribution"
date:   2020-12-12 14:14:14 +0100
tags: frequentist statistics exponential distribution
categories: frequentist statistics
math: true
---

In this post we show two different exact confidence intervals for the rate parameter $$\lambda$$ of exponentially-distributed random variables. One useful, and one less so.
First, suppose we have $$X_1, \ldots, X_n \sim \text{Exp}(\lambda)$$, with $$\lambda$$ unknown. Based on our sample, we want to get a confidence interval for $$\lambda$$ of confidence level $$1-\alpha$$.
To create an exact confidence interval, we need a pivot: a function $$T(X_1, \ldots, X_n, \theta)$$ of the observations and the parameter, whose probability does not depend on $$\theta$$. Rather, this function follows some known distribution from which we can calculate the quantiles. The examples below will show you how such a pivot can lead to a confidence interval.

We start by looking at the distribution function of the exponential distribution for some $$X \sim \text{Exp}(\lambda)$$:

$$ 
F(x;\lambda) = 1 - e^{-\lambda x}, \qquad x \geq 0. 
$$

If we are smart we can spot that if we fill in $$x/\lambda$$ instead of $$x$$, we end up with a standard Exponential distribution. This is useful, because it shows us how to transform our samples such that they become standard distributed:

$$
1 - e^{-x} = P(X \leq x/\lambda) = P(\lambda X \leq x) \quad \Rightarrow \quad \lambda X \sim \text{Exp}(1).
$$

Now often to create a confidence interval we look at what the distribution is of the mean of the data $$\overline{X}$$ of the sum of the data $$\sum_{i=1}^n X_i$$. 

The first confidence interval is the most convenient one as it makes use of the chi-squared distribution, for which cumulative probability tables exist in books and online. Note that if we take $$\lambda = \frac{1}{2}$$ then $$X \sim \chi_2^2$$ which we can see in the following way:

$$ 
f_X(x/2) = \frac{e^{-\frac{x}{2}}}{2} = \frac{x^{\frac{2}{2}-1}e^{-\frac{x}{2}}}{2^{\frac{2}{2}}\Gamma(\frac{2}{2})}
$$

where $$f_X$$ is the density function of $$X$$ and on the right side we have the density function of the Gamma distribution with $$k = 2$$.
In the same way as above, we note that $$2\lambda X_i \sim \text{Exp}(1/2)$$ for each $$X_i$$. We can also make use of the sum, since the sum of $n$ chi-squared distributions with $k$ degrees of freedom each, is again a chi-squared distribution with $$nk$$ degrees of freedom. This can be seen by writing each chi-squared distribution as the sum of standard normal random variables, and taking the sum. So we end up with 

$$
2\lambda X_1 + \ldots + 2\lambda X_n = 2\lambda \sum_{i=1}^n X_i \sim \chi_{2n}^2.
$$

This allows us to create a confidence interval:

$$
\begin{align*} 
1 - \alpha &= P_\lambda\Bigg[\chi^2_{2n,\alpha/2} \leq 2\lambda \sum_{i=1}^n X_i \leq \chi^2_{2n, 1-\alpha/2}\Bigg] \\ &= P_\lambda\Bigg[\frac{\chi^2_{2n,\alpha/2}}{2\sum_{i=1}^n X_i} \leq \lambda  \leq \frac{\chi^2_{2n, 1-\alpha/2}}{2\sum_{i=1}^n X_i}\Bigg]
\end{align*}
$$

So that we have the following confidence interval $$[\frac{\chi^2_{2n,\alpha/2}}{2\sum_{i=1}^n X_i}, \frac{\chi^2_{2n, 1-\alpha/2}}{2\sum_{i=1}^n X_i}]$$, and we are done. 

For the second confidence interval we look at the order statistic $$X_{(1)}$$ of our sample, the smallest value we observed. Let's find its distribution:

$$ 
\begin{align*}
P(X_{(1)} \leq x) &= 1 - P(X_1, \ldots, X_n > 1) \\ &= 1 - \prod_{i=1}^n P(X_i > 1) \\ &= 1 - \prod_{i=1}^n (1 - P(X_i \leq 1)) \\ &= 1 - \prod_{i=1}^n e^{-\lambda x} \\ &= 1 - e^{-n\lambda x}.
\end{align*}
$$

which means $$X_{(1)}$$ is again exponentially distributed with parameter $$n\lambda$$! Again, we can take $$n\lambda X_{(1)}$$, which then follows a standard exponential distribution. Now if we let $$E_{1-\alpha}$$ be the value such that $$P(X \leq E_{\alpha}) = 1-\alpha$$, then we can create the following confidence interval: 

$$
P(n\lambda X_{(1)} \leq E_{1-\alpha}) = P(\lambda \leq E_{1-\alpha}/(nX_{(1)}) = 1-\alpha.
$$

Note that while correct, this confidence interval is way less useful than our first one, which we can see by performing some simulations. We can also start it from some value in $$[0, \alpha]$$ so that we get $$P(E_{\beta}/(nX_{(1)}) \leq \lambda \leq E_{1-\gamma}/(nX_{(1)})) = 1-\alpha$$ as long as $$\beta + \gamma = \alpha$$.

{% highlight R %}
lambda = 12
N = 50 # amount of samples
Q = 10000 # amount of experiments
alpha = 0.05

count_in_ci = rep(0, Q)
q1 = qchisq(alpha/2, 2*N)
q2 = qchisq(1-alpha/2, 2*N)
for (i in 1:Q) { 
	sample = rexp(N, lambda)
	lower_ci = q1/(2*sum(sample))
	upper_ci = q2/(2*sum(sample))
	count_in_ci[i] = (lower_ci <= lambda) & (lambda <= upper_ci)
}
lower_ci # [1] 7.986
upper_ci # [1] 13.941
sum(count_in_ci)/Q # [1] 0.947

count_in_ci = rep(0, Q)
q = qexp(1-alpha)
for (i in 1:Q) { 
	sample = rexp(N, lambda)
	smallest = min(sample)
	upper_ci = q/(N*smallest)
	count_in_ci[i] = (lambda <= upper_ci)
}
upper_ci # [1] 24.960
sum(count_in_ci)/Q # [1] 0.953

{% endhighlight %}

Which shows us that both confidence intervals indeed perform as expected.