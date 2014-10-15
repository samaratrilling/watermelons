Project 3: Watermelon

The genetics of watermelon breeding are quite interesting. One class of natural watermelons is diploid, meaning that (like humans) there are two copies of each chromosome, one inherited from each parent. A second class is tetraploid: two copies of each chromosome are inherited from each parent, making four copies in total.

An intriguing result of breeding a diploid plant with a tetraploid plant is that the progeny is triploid. A side-effect of triploidy is that the plant produces almost no seeds. Interbreeding would be selected against in the wild, since seedless watermelons would not reproduce. Nevertheless, humans can overcome this barrier by explicitly creating conditions in which diploid and tetraploid plants pollinate each other. In fact, humans are motivated to do so because seedless watermelons are just so much easier to eat!

Your job is to plant a field with diploid or tetraploid seeds in such a way as to maximize the output of fruit from a field. The metric rewards the production of seedless watermelons: The reward for producing a seedless watermelon is one dollar, whereas the reward for producing a watermelon containing seeds is s dollars, where 0<s<1. s is a parameter that we will vary.

You will be given a rectangular field of dimensions L by W meters. This field has a number of pre-existing trees within it, whose coordinates are specified. A seed must be planted in such a way that it has exclusive access to soil within one meter. That means each seed must be at least one meter from any field boundary, and at least two meters from all other seeds and trees. You specify the location of each seed, and its ploidy.

Bees pollinate the crop by flying between plants. Since bees are more likely to fly between nearby plants, the probability of pollination is higher when plants are close together. We model this process as an inverse-square law, so that the probability of a plant P being pollinated by a plant Q is 

\begin{displaymath}
\frac{\vert PQ\vert^{-2}}{\sum_{R \not = P} \vert PR\vert^{-2}}
\end{displaymath}

where |PQ| is the Euclidean distance from P to Q, and R ranges over all plants in the field other than P.
For each plant P, we calculate the probability XP that it will be pollinated by a plant of the opposite ploidy. The total score in dollars is then 

\begin{displaymath}
\sum_P X_P + s(1-X_P)
\end{displaymath}

This function rewards both packing plants close together (so that many fit within the field) as well as arranging the ploidies to promote cross-pollination. In some cases, these submetrics conflict (when?) and you will need to make an informed decision about the trade-offs.
We will supply a simulator that will check all of the constraints mentioned above, and calculate the overall score of your configuration. At the end of the class we'll run a tournament using various parameter settings and tree layouts, including some that will be unseen.

Some initial things to think about:

Try thinking about simpler configurations first, such as those without trees.
A one-dimensional version of the problem can be specified by setting W=2. Does this version of the problem have an exact solution?
Packing equal circles into squares has been previously studied. See this reference for example. How might you use such information?
