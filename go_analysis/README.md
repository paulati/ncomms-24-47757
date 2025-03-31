## Go term analysis

The binomial distribution models the probabilities of having a certain number of successes among N identical trials (each having p as the probability of success). It has discrete values, counts the number of successes in yes/no-type experiments. Each of these experiments, also called Bernoulli trials, results in either success or failure.

This distribution is the one that models the region-centric problem solved by rGREAT.

rGREAT null assumption is 'genomic regions are uniformly distributed on the genome', and the enrichment test is applied as follows. For a specific biological term represented as a gene set, denote the fraction of its associated functional domains in the genome as p, the total number of input regions as N, the observed number of input regions that fall in the associated domains as n and the corresponding random variable as X, then X follows a binomial distribution: X ∼ B(p, N) and the p-value of the enrichment is calculated as Pr(X ≥ n). In our problem, X is a random variable counting the overlap between (accelerated or conserved) genomic regions and a specific go-term associated functional domains (observed regions hits), n is the actual value of X, and N is the accelerated elements count.

We want to evaluate go term enrichment for the set of regions associated with accelerated elements. We have an extra issue trying to know if the significance is associated with acceleration or if it is just a consequence of enrichment in the conserved regions set. (Remember that accelerated regions are a subset of conserved regions.)

The strategy to answer this question is to compare results obtained for the set of accelerated regions with those obtained from equivalent sets of conserved regions. We sampled conserved regions, building 5,000 sets with the same size and chromosome distribution as the accelerated regions set, and evaluated rGREAT on each of these sets.

From RGREAT source code available at <https://github.com/jokergoo/rGREAT/blob/master/R/great_local.R>

We have:

<https://github.com/jokergoo/rGREAT/blob/master/R/great_local.R#L608>

`prop = width_fgr/background_total_length`

<https://github.com/jokergoo/rGREAT/blob/master/R/great_local.R#L618>

`gr = intersect(gr, background)`

`background` is the whole genome, then `gr` are the input regions (accelerated/conserved)

<https://github.com/jokergoo/rGREAT/blob/master/R/great_local.R#L523>

`n_total = length(gr)`

`n_total` is the same as N in the previous paragraph, the count of input elements.

<https://github.com/jokergoo/rGREAT/blob/master/R/great_local.R#L618>

`p = pbinom(q = (n_hits - 1), size = n_total, prob = prop, log.p = TRUE, lower.tail = FALSE)`

<https://github.com/jokergoo/rGREAT/blob/master/R/great_local.R#L643>

`p_value = p, p_adjust = p.adjust(p, "BH"),`

We observe that `prop` value is approximately the same for both a conserved region subset and the accelerated regions set. By design, `width_fgr` values are also similar between the two sets, and the `background_total_length` is exactly the same (whole genome)

`gr` is the input regions set (conserved/accelerated), and the intersection of regions with the whole genome returns the same gr value. `n_total` value is the same in both sets by design. `p` depends on `n_total`, `prop`, and `n_hits`. Since `n_total` and `prop` values are the same for the accelerated regions set and the conserved regions subsets, in these specific conditions p depends only on `n_hits`.

Based on these considerations, we propose the following strategy to select go terms associated with accelerated regions. We build an empirical distribution of `n_hits` values for a specific go term from rGREAT results of each of the 5,000 conserved regions subsets. Then we calculate `n_hits` value from rGREAT result of the accelerated regions set (acc_observed_region_hits). We calculate the proportion of random n_hits observed that are greater than acc_observed_region_hits. If this proportion is less than 0.05, then this set of accelerated regions differs from the random sets of conserved regions evaluating this specific go term and we can assume that the behavior of the accelerated regions set is independent of its conserved nature for this term. We report a go terms as associated with the accelerated regions set if the behavior of the accelerated regions set is independent from its conserved nature and adjusted p-value that is the result of rGREAT on accelerated regions set is less than 0.05.

References:

<https://academic.oup.com/bioinformatics/article/39/1/btac745/6832038>

<https://github.com/jokergoo/rGREAT>
