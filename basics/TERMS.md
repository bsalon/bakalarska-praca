## Cryptography terminology
  
  
Good site for [basic terminology](http://www.aspencrypt.com/crypto101_terminology.html).
&nbsp;
  
* cryptographic primitives

  are well-established, low-level cryptographic algorithms that are frequently used to build\
  cryptographic protocols for computer security systems.
  * Examples: one-way hash function, symmetric key criptography, public-key crypthography, digital\
    signatures, mix network, pribate information retrieval, commitment scheme, cryptographically secure\
    pseudorandom number generator\
&nbsp;
&nbsp;
* **distinguisher**

  is an arbitrary algorithm. In fact, we do NOT want to formalize anything about the distinguisher\
  (except that its output is a single bit, although we don't even really need to do this). In\
  definitions, we require that no distinguisher should succeed with non-negligible probability. So,\
  this should hold for any algorithm.\
  \
  think of a distinguisher as an adversarially-chosen statistical experiment that attempts to support\
  or refute some hypothesis, that again is adversarially selected. This means that the honest party\
  doesn't get to assume what meaning the adversary assigns to 0 or 1; their PRG has to behave like\
  a random distribution no matter what the adversary does (with the proviso that the adversary's\
  experiment must complete within polynomial time).\
&nbsp;
&nbsp;
* block cipher

  is a deterministic algorithm operating on fixed-length groups of bits, called blocks. It uses an\
  unvarying transformation, that is, it uses a symmetric key. They are specified elementary components\
  in the design of many cryptographic protocols and are widely used to implement the encryption of\
  large amounts of data, including data exchange protocols.\
&nbsp;
&nbsp;
* stream ciphers

  is a symmetric key cipher where plaintext digits are combined with a pseudorandom cipher digit stream\
  (keystream). In a stream cipher, each plaintext digit is encrypted one at a time with the\
  corresponding digit of the keystream, to give a digit of the ciphertext stream. Since encryption of\
  each digit is dependent on the current state of the cipher, it is also known as state cipher. In\
  practice, a digit is typically a bit and the combining operation is an exclusive-or (XOR).\
&nbsp;
&nbsp;
* cryptoanalysis

  is the study of analyzing information systems in order to study the hidden aspects of the systems. It\
  is used to breach cryptographic security systems and gain access to the contents of encrypted\
  messages, even if the cryptographic key is unknown. In addition to mathematical analysis of\
  cryptographic algorithms, cryptanalysis includes the study of side-channel attacks that do not target\
  weaknesses in the cryptographic algorithms themselves, but instead exploit weaknesses in their\
  implementation.\
&nbsp;
&nbsp;
* statistical test

  provides a mechanism for making quantitative decisions about a process or processes. The intent is to\
  determine whether there is enough evidence to "reject" a conjecture or hypothesis about the process.\
&nbsp;
&nbsp;
* Monobit test

  treats each output bit of the random number generator as a coin flip test, and determine if the\
  observed number of heads and tails are close to the expected 50% frequency. The number of heads in\
  a coin flip trail forms a binomial distribution. It tests statistical randomness.\
&nbsp;
&nbsp;
* Fourier Transform

  is a mathematical transform that decomposes a function (often a function of time, or a signal) into\
  its constituent frequencies, such as the expression of a musical chord in terms of the volumes and\
  frequencies of its constituent notes. The term Fourier transform refers to both the frequency domain\
  representation and the mathematical operation that associates the frequency domain representation\
  to a function of time.\
&nbsp;
&nbsp;
* hamming weight

  of a string is the number of symbols that are different from the zero-symbol of the alphabet used. It\
  is thus equivalent to the Hamming distance from the all-zero string of the same length. For the most\
  typical case, a string of bits, this is the number of 1's in the string, or the digit sum of the\
  binary representation of a given number and the ℓ₁ norm of a bit vector. In this binary case, it is\
  also called the population count, popcount, sideways sum, or bit summation.\
&nbsp;
&nbsp;
* **plaintext - primitive input**

  usually means unencrypted information pending input into cryptographic algorithms, usually encryption\
  algorithms. Cleartext usually refers to data that is transmitted or stored unencrypted ("in clear").\
&nbsp;
&nbsp;
* **ciphertext - primitive output**

  is the result of encryption performed on plaintext using an algorithm, called a cipher. Ciphertext is\
  also known as encrypted or encoded information because it contains a form of the original plaintext\
  that is unreadable by a human or computer without the proper cipher to decrypt it. Decryption, the\
  inverse of encryption, is the process of turning ciphertext into readable plaintext.\
&nbsp;
&nbsp;
* **key (bits)**

  is a piece of information (a parameter) that determines the functional output of a cryptographic\
  algorithm. For encryption algorithms, a key specifies the transformation of plaintext into\
  ciphertext, and vice versa for decryption algorithms. Keys also specify transformations in other\
  cryptographic algorithms, such as digital signature schemes and message authentication codes.\
&nbsp;
&nbsp;
* linear cryptoanalysis

  is a general form of cryptanalysis based on finding affine approximations to the action of a cipher.\
  Attacks have been developed for block ciphers and stream ciphers. Linear cryptanalysis is one of the\
  two most widely used attacks on block ciphers; the other being differential cryptanalysis.\
&nbsp;
&nbsp;
* differential cryptoanalysis

  is a general form of cryptanalysis applicable primarily to block ciphers, but also to stream ciphers\
  and cryptographic hash functions. In the broadest sense, it is the study of how differences in\
  information input can affect the resultant difference at the output. In the case of a block cipher,\
  it refers to a set of techniques for tracing differences through the network of transformation,\
  discovering where the cipher exhibits non-random behavior, and exploiting such properties to recover\
  the secret key (cryptography key).\
&nbsp;
&nbsp;
* algebraic cryptoanalysis

  is the process of breaking codes by solving polynomialsystems of equations\
&nbsp;
&nbsp;
* **Z-score [standard score]**

  is the number of standard deviations by which the value of a raw score (i.e., an observed value or\
  data point) is above or below the mean value of what is being observed or measured. Raw scores above\
  the mean have positive standard scores, while those below the mean have negative standard scores. It\
  is calculated by subtracting the population mean from an individual raw score and then dividing the\
  difference by the population standard deviation. This process of converting a raw score into\
  a standard score is called standardizing or normalizing (however, "normalizing" can refer to many\
  types of ratios).\
&nbsp;
  The Z-score represents the distance from the mean μ in units of σ. The binomial distribution B(n,p)\
  is approximated by N(μ,σ²), with the parameters μ=np and σ²=np(1−p) i.e.Z-score of binomially\
  (B(n,p)) distributed #1 is computed as Z-score = (#1 − pn) / √(p(1 − p)n) = #1 − μσ\
&nbsp;
&nbsp;
* **boolean function**

  is a function whose arguments, as well as the function itself, assume values from a two-element set\
  (usually {0,1}). As a result, it is sometimes referred to as a "switching function". A Boolean\
  function takes the form f : {0, 1}ⁿ → {0, 1}, where {0, 1} is called a Boolean domain and n is a\
  non-negative integer called the arity of the function. In the case where n = 0 the "function"\
  is essentially a constant element of {0, 1}.\Every n-ary Boolean function can be expressed as a\
  propositional formula in n variables x_1, ..., x_n, and two propositional formulas are logically\
  equivalent if and only if they express the same Boolean function. There are 2^(2ⁿ) n-ary functions\
  for every n.
&nbsp;
&nbsp;
* histogram

  is an approximate representation of the distribution of numerical data.\
&nbsp;
&nbsp;
* **null hypothesis**

  is a general statement or default position that there is no relationship between two measured\
  phenomena or no association among groups. Testing (rejecting or not rejecting) the null hypothesis\
  and thus concluding that there are (or there are not) grounds for believing that there is\
  a relationship between two phenomena (e.g., that a potential treatment has a measurable effect) is\
  a central task in the modern practice of science; the field of statistics, more specifically\
  hypothesis testing, gives precise criteria for rejecting or not rejecting a null hypothesis within\
  a confidence level. In BoolTest - "data is random".\
&nbsp;
&nbsp;
* **observed test statistic s_obs**

  is a statistic (a quantity derived from the sample) used in statistical hypothesis testing.\
  A hypothesis test is typically specified in terms of a test statistic, considered as a numerical\
  summary of a data-set that reduces the data to one value that can be used to perform the hypothesis\
  test. In general, a test statistic is selected or defined in such a way as to quantify, within\
  observed data, behaviours that would distinguish the null from the alternative hypothesis, where such\
  an alternative is prescribed, or that would characterize the null hypothesis if there is no\
  explicitly stated alternative hypothesis. An important property of a test statistic is that its\
  sampling distribution under the null hypothesis must be calculable, either exactly or approximately,\
  which allows p-values to be calculated.\
&nbsp;
&nbsp;
* **p-value**

  In null hypothesis significance testing, the p-value is the probability of obtaining test results\
  at least as extreme as the results actually observed, under the assumption that the null hypothesis\
  is correct.(In the case of a composite null hypothesis, the largest such probability allowed under\
  the null hypothesis is taken.) A very small p-value means that such an extreme observed outcome\
  would be very unlikely under the null hypothesis. Reporting p-values of statistical tests is common\
  practice in academic publications of many quantitative fields. Since the precise meaning of p-value\
  is hard to grasp, misuse is widespread and has been a major topic in metascience.\
&nbsp;
  The p-value represents the probability that more extreme results are obtained (for the true\
  hypothesis) than we observed (s_obs). In our case, p-value represents the probability that a perfect\
  random number generator would produce less random sequences than the sequence being tested. The\
  p-value is computed from the observed test statistic s_obs and the reference distribution D or its\
  close approximation. The null distribution of many tests is binomial distribution B(n,p). It is\
  approximated well (for n>10p and n·(1−p) > 10) by normal distribution N(μ,σ²). Normal distribution is\
  symmetric around mean μ and therefore p-value is computed as an area under bell curve in both tails.\
  Sometimes Z-score is computed instead of p-value since they are related:\
    p-value = erfc(Z-score√2)\
  and computation of Z-score is simpler and faster.\
&nbsp;
  The p-value represents the probability that RNG would generate data with more extreme test statistic\
  than s_obs. A p-value can be computed as an area below the normal distribution in the tail bounded by\
  the observed test statistic s_obs.
&nbsp;
&nbsp;
* Test-U01 (state-of-art batt.)

  is a software library, implemented in the ANSI C language, that offers a collection of utilities for\
  the empirical randomness testing of random number generators (RNGs).\
&nbsp;
&nbsp;
* NIST STS

  A Statistical Test Suite for Random and Pseudorandom Number Generators for Cryptographic\
  Applications.\
&nbsp;
&nbsp;
* Dieharder

  are a battery of statistical tests for measuring the quality of a random number.
&nbsp;
&nbsp;
* profiling

  is a form of dynamic program analysis that measures, for example, the space (memory) or time\
  complexity of a program, the usage of particular instructions, or the frequency and duration of\
  function calls. Most commonly, profiling information serves to aid program optimization.\
&nbsp;
&nbsp;
* Type I and Type II errors

  In statistical hypothesis testing, a type I error is the rejection of a true null hypothesis (also\
  known as a "false positive" finding or conclusion; example: "an innocent person is convicted"),\
  while a type II error is the non-rejection of a false null hypothesis (also known as a\
  "false negative" finding or conclusion; example: "a guilty person is not convicted"). Much of\
  statistical theory revolves around the minimization of one or both of these errors, though the\
  complete elimination of either is a statistical impossibility for non-deterministic algorithms. By\
  selecting a low threshold (cut-off) value and modifying the alpha (p) level, the quality of the\
  hypothesis test can be increased. Intuitively, type I errors can be thought of as errors of\
  commission, and type II errors as errors of omission.
  * Null hypothesis - "data is random"\
    Type I error    - truly random data rejected\
    Type II error   - non-random data not rejected\
&nbsp;
&nbsp;
* standard security parameters

  is a way of measuring of how "hard" it is for an adversary to break a cryptographic scheme. There are\
  two main types of security parameter: computational and statistical, often denoted by κ and λ,\
  respectively. Roughly speaking, the computational security parameter is a measure for the input size\
  of the computational problem on which the cryptographic scheme is based, which determines its\
  computational complexity, whereas the statistical security parameter is a measure of the probability\
  with which an adversary can break the scheme (whatever that means for the protocol).\
&nbsp;
&nbsp;
* distiguishing attack

  is any form of cryptanalysis on data encrypted by a cipher that allows an attacker to distinguish the\
  encrypted data from random data. Modern symmetric-key ciphers are specifically designed to be immune\
  to such an attack. In other words, modern encryption schemes are pseudorandom permutations and are\
  designed to have ciphertext indistinguishability. If an algorithm is found that can distinguish the\
  output from random faster than a brute force search, then that is considered a break of the cipher.\
  A similar concept is the known-key distinguishing attack, whereby an attacker knows the key and can\
  find a structural property in cipher, where the transformation from plaintext to ciphertext is\
  not random.\
&nbsp;
&nbsp;
* avalanche effect

  is the desirable property of cryptographic algorithms, typically block ciphers and cryptographic hash\
  functions, wherein if an input is changed slightly (for example, flipping a single bit), the output\
  changes significantly (e.g., half the output bits flip). In the case of high-quality block ciphers,\
  such a small change in either the key or the plaintext should cause a drastic change in the\
  ciphertext.\
&nbsp;
&nbsp;
* strict avalanche critetion

  is a formalization of the avalanche effect. It is satisfied if, whenever a single input bit is\
  complemented, each of the output bits changes with a 50% probability.\
&nbsp;
&nbsp;
* rank of a matrix

  is the dimension of the vector space generated (or spanned) by its columns. This corresponds to the\
  maximal number of linearly independent columns of A. This, in turn, is identical to the dimension of\
  the vector space spanned by its rows. Rank is thus a measure of the "nondegenerateness" of the system\
  of linear equations and linear transformation encoded by A.\
&nbsp;
&nbsp;
* frequency test

  is used to test the randomness of a sequence of zeroes and ones. The test is based on the proportion\
  of zeroes and ones. Specifically, it tests the closeness of the proportion of ones to 0.5. The\
  frequency within a block test is a refinement that tests the proportion of ones within M-value blocks.\
&nbsp;
&nbsp;
* **normal distribution**

  is a type of continuous probability distribution for a real-valued random variable. The general form\
  of its probability density function is f(x) = 1 / (o sqrt(2pi) ) * e ^ (-1 / 2 * ((x - u) / o)²). The\
  parameter μ is the mean or expectation of the distribution (and also its median and mode), while the\
  parameter σ is its standard deviation. The variance of the distribution is σ². A random variable with\
  a Gaussian distribution is said to be normally distributed, and is called a normal deviate.\
&nbsp;
&nbsp;
* null distribution

  is the probability distribution of the test statistic when the null hypothesis is true. For example,\
  in an F-test, the null distribution is an F-distribution. Null distribution is a tool scientists\
  often use when conducting experiments. The null distribution is the distribution of two sets of data\
  under a null hypothesis. If the results of the two sets of data are not outside the parameters of the\
  expected results, then the null hypothesis is said to be true.\
&nbsp;
&nbsp;
* mean

  For a data set, the arithmetic mean, also called the expected value or average, is the central value\
  of a discrete set of numbers: specifically, the sum of the values divided by the number of values.\
  The arithmetic mean of a set of numbers x1, x2, ..., xn is typically denoted by x¯.\
&nbsp;
&nbsp;
* standard deviation

  is a measure of the amount of variation or dispersion of a set of values. A low standard deviation\
  indicates that the values tend to be close to the mean (also called the expected value) of the set,\
  while a high standard deviation indicates that the values are spread out over a wider range.\
&nbsp;
&nbsp;
* error function

  is a complex function of a complex variable defined as: erf z = ( 2 / sqrt(π) ) ∫<0;z> e^(−t^2) dt.\
  This integral is a special (non-elementary) and sigmoid function that occurs often in probability,\
  statistics, and partial differential equations. In many of these applications, the function argument\
  is a real number. If the function argument is real, then the function value is also real. When the\
  results of a series of measurements are described by a normal distribution with standard deviation σ\
  and expected value 0, then erf ( a / σ sqrt(2) ) is the probability that the error of a single\
  measurement lies between −a and +a, for positive a. This is useful, for example, in determining the\
  bit error rate of a digital communication system.\
  Two closely related functions are the complementary error function (erfc) defined as:\
   erfc z = 1 − erf z\
  and the imaginary error function (erfi) defined as:\
   erfi z =  i er(iz)\
  where i is the imaginary unit.\
&nbsp;
&nbsp;
* **binomical distribution**

  with parameters n and p, denoted Bin (n ,p) is the discrete probability distribution of the number of\
  successes in a sequence of n independent experiments, each asking a yes–no question, and each with\
  its own boolean-valued outcome: success/yes/true/one (with probability p), or failure/no/false/zero\
  (with probability q = 1 − p). A single success/failure experiment is also called a Bernoulli trial\
  or Bernoulli experiment, and a sequence of outcomes is called a Bernoulli process; for a single trial\
  (i.e., n = 1), the binomial distribution is a Bernoulli distribution. The binomial distribution is\
  the basis for the popular binomial test of statistical significance.\
&nbsp;
&nbsp;
* distributions

  are objects that generalize the classical notion of functions in mathematical analysis. Distributions\
  make it possible to differentiate functions whose derivatives do not exist in the classical sense. In\
  particular, any locally integrable function has a distributional derivative.\
&nbsp;
&nbsp;
* **maximal Z-SCORE**

  In order to interpret result (Z-SCORE) of our tool we have to find expected distribution f_max of\
  Z-SCORE for the null hypothesis i.e. random data. The value of Z-SCORE is determined by the boolean\
  functions constructed in two phases using setting deg,m,k,t. The value Z-SCORE is computed as\
  a maximal Z-score for a set of Z-scores corresponding to constructed boolean functions. The\
  theoretical assessment of the distribution f_max should follow the process of the construction of\
  the best distinguisher with maximal Z-score (in absolute value). It is hard to theoretically derive\
  probability density function (pdf) f_max for settings with k>1. But for the settings with k=1\
  the pdf f_max can be derived much easier.\
&nbsp;
&nbsp;
* normalization (statistics)

  can have a range of meanings. In the simplest cases, normalization of ratings means adjusting values\
  measured on different scales to a notionally common scale, often prior to averaging. In more\
  complicated cases, normalization may refer to more sophisticated adjustments where the intention is\
  to bring the entire probability distributions of adjusted values into alignment.\
&nbsp;
&nbsp;
* term

  is either a single number or variable, or the product of several numbers or variables. Terms are\
  separated by a + or - sign in an overall expression.\
&nbsp;
&nbsp;
* monomial

  is, roughly speaking, a polynomial which has only one term.\
  x1 XOR x2 XOR x3 has three monomials.\
&nbsp;
&nbsp;
* basis

  In mathematics, a set B of elements (vectors) in a vector space V is called a basis, if every element\
  of V may be written in a unique way as a (finite) linear combination of elements of B. The\
  coefficients of this linear combination are referred to as components or coordinates on B of\
  the vector. The elements of a basis are called basis vectors. Equivalently B is a basis if its\
  elements are linearly independent and every element of V is a linear combination of elements of B. In\
  more general terms, a basis is a linearly independent spanning set.\
&nbsp;
&nbsp;
* algebraic normal form

  is a way of writing logical formulas in one of three subforms:\
   The entire formula is purely true or false:\
      1\
      0\
   One or more variables are ANDed together into a term, then one or more terms are XORed together\
   into ANF. No NOTs are permitted:\
      a ⊕ b ⊕ ab ⊕ abc\
    or in standard propositional logic symbols:\
      a ⊻ b ⊻ ( a ∧ b ) ⊻ ( a ∧ b ∧ c )\
   The previous subform with a purely true term:\
      1 ⊕ a ⊕ b ⊕ ab ⊕ abc\
&nbsp;
&nbsp;
* entropy

  of a random variable is the average level of "information", "surprise", or "uncertainty" inherent in
  the variable's possible outcomes. As an example, consider a biased coin with probability p of landing\
  on heads and probability 1-p of landing on tails. The maximum surprise is for p = 1/2, when there is\
  no reason to expect one outcome over another, and in this case a coin flip has an entropy of one bit.\
  The minimum surprise is when p = 0 or p = 1, when the event is known and the entropy is zero bits.\
  Other values of p give different entropies between zero and one bits.\
&nbsp;
&nbsp;
* probability density function

  of a continuous random variable, is a function whose value at any given sample (or point) in the\
  sample space (the set of possible values taken by the random variable) can be interpreted as\
  providing a relative likelihood that the value of the random variable would equal that sample. In\
  other words, while the absolute likelihood for a continuous random variable to take on any particular\
  value is 0 (since there are an infinite set of possible values to begin with), the value of the PDF\
  at two different samples can be used to infer, in any particular draw of the random variable, how\
  much more likely it is that the random variable would equal one sample compared to the other sample.\
  In a more precise sense, the PDF is used to specify the probability of the random variable falling\
  within a particular range of values, as opposed to taking on any one value.\
&nbsp;
&nbsp;
* Galois field [GF(2) -> x^2 + x + 1]

  is a field that contains a finite number of elements. As with any field, a finite field is a set on\
  which the operations of multiplication, addition, subtraction and division are defined and satisfy\
  certain basic rules. The most common examples of finite fields are given by the integers mod p when p\
  is a prime number.\
&nbsp;
&nbsp;
* tweakable polynomials over GF(2)

  Almost any cryptographic scheme can be described by tweakable polynomials over GF(2), which contain\
  both secret variables (e.g., key bits) and public variables (e.g., plaintext bitsor IV bits). The\
  cryptanalyst is allowed to tweak the polynomials by choosing arbitrary values for the public\
  variables, and his goal is to solve the resultant system of polynomial equations in terms of their\
  common secret variables. In this paper we develop a new technique (called a cube attack) for solving\
  such tweakable polynomials, which is a major improvement over several previously published attacks\
  of the same type.\
&nbsp;
&nbsp;
* cube attack

  Over a general field GF(p^k) with p > 2, the correct way to apply cube attacks is to alternately add\
  and subtract the outputs of the master polynomial with public inputs that range only over the two\
  values 0 and 1 (and not over all their possible values of 0,1,2, ..., p−1), where the sign is\
  determinedby the sum (modulo 2) of the vector of assigned values.\
&nbsp;
&nbsp;
* IV (initialization vector)

  is a fixed-size input to a cryptographic primitive that is typically required to be random or\
  pseudorandom. Randomization is crucial for encryption schemes to achieve semantic security, a property\
  whereby repeated usage of the scheme under the same key does not allow an attacker to infer\
  relationships between segments of the encrypted message. For block ciphers, the use of an IV\
  is described by the modes of operation. Randomization is also required for other primitives, such as\
  universal hash functions and message authentication codes based thereon.\
&nbsp;
&nbsp;
* linearization technique

  is finding the linear approximation to a function at a given point. The linear approximation of a\
  function is the first order Taylor expansion around the point of interest. In the study of dynamical\
  systems, linearization is a method for assessing the local stability of an equilibrium point of\
\  a system of nonlinear differential equations or discrete dynamical systems.\
&nbsp;
* superset

  A set A is a superset of another set B if all elements of the set B are elements of the set A. The\
  superset relationship is denoted as A ⊃ B.\
&nbsp;
&nbsp;
* index set

  is a set whose members label (or index) members of another set. For instance, if the elements of\
  a set A may be indexed or labeled by means of the elements of a set J, then J is an index set. The\
  indexing consists of a surjective function from J onto A and the indexed collection is typically\
  called an (indexed) family, often written as {A_j} j ∈ J.\
&nbsp;
&nbsp;
* maxterm

  is a sum that is a logical OR of a set of variables where each individual variable only appears once\
  in the sum, either in complemented or uncomplemented form, so that the value of the sum becomes 0.\
  All the maxterms in a product of maxterms should have the same variables, although each maxterm\
  should differ from every other one by the pattern of complementation of those variables.\
&nbsp;
&nbsp;
* FFT (fast fourier transform)

  is an algorithm that computes the discrete Fourier transform (DFT) of a sequence, or its inverse\
  (IDFT). Fourier analysis converts a signal from its original domain (often time or space) to a\
  representation in the frequency domain and vice versa. The DFT is obtained by decomposing a sequence\
  of values into components of different frequencies. This operation is useful in many fields, but\
  computing it directly from the definition is often too slow to be practical. An FFT rapidly computes\
  such transformations by factorizing the DFT matrix into a product of sparse (mostly zero) factors.\
&nbsp;
&nbsp;
* free term

  A term is said to be variable-free if no variables occur in it.\
&nbsp;
&nbsp;
* nonsingular(invertible) matrix

  there exists an n-by-n square matrix B such that, AB = BA = I_n, where I_n denotes the n-by-n\
  identity matrix and the multiplication used is ordinary matrix multiplication. If this is the case,\
  then the matrix B is uniquely determined by A, and is called the (multiplicative) inverse of A,\
  denoted by A^-1.
&nbsp;
&nbsp;
* singular matrix

  is singular if and only if its determinant is zero. Singular matrices are rare in the sense that if\
  a square matrix's entries are randomly selected from any finite region on the number line or complex\
  plane, the probability that the matrix is singular is 0, that is, it will "almost never" be singular.\
&nbsp;
&nbsp;
* PCP theorem

  states that every decision problem in the NP complexity class has probabilistically checkable proofs\
  (proofs that can be checked by a randomized algorithm) of constant query complexity and logarithmic\
  randomness complexity (uses a logarithmic number of random bits). The PCP theorem says that for some\
  universal constant K, for every n, any mathematical proof of length n can be rewritten as a different\
  proof of length poly(n) that is formally verifiable with 99% accuracy by a randomized algorithm that\
  inspects only K letters of that proof.\
&nbsp;
&nbsp;
* CLR test

  conditional likelihood ratio test\
&nbsp;
&nbsp;
* message authentication code

  is a short piece of information used to authenticate a message—in other words, to confirm that the\
  message came from the stated sender (its authenticity) and has not been changed. The MAC value\
  protects both a message's data integrity as well as its authenticity, by allowing verifiers (who also\
  possess the secret key) to detect any changes to the message content.\
&nbsp;
&nbsp;
* Chosen-IV (attack)

  Stream ciphers combine a secret key with an agreed initialization vector (IV) to produce a\
  pseudo-random sequence which from time-to-time is re-synchronized. A "Chosen IV" attack relies on\
  finding particular IV's which taken together probably will reveal information about the secret key.\
  Typically multiple pairs of IV are chosen and differences in generated key-streams are then analysed\
  statistically for a linear correlation and/or an algebraic boolean relation (see also Differential\
  cryptanalysis). If choosing particular values of the initialization vector does expose a non-random\
  pattern in the generated sequence, then this attack computes some bits and thus shortens the\
  effective key length. A symptom of the attack would be frequent re-synchronisation. Modern stream\
  ciphers include steps to adequately mix the secret key with an initialization vector, usually by\
  performing many initial rounds.\
&nbsp;
&nbsp;
* linear feedback shift register

  is a shift register whose input bit is a linear function of its previous state. The most commonly\
  used linear function of single bits is exclusive-or (XOR). Thus, an LFSR is most often a shift\
  register whose input bit is driven by the XOR of some bits of the overall shift register value.\
&nbsp;
&nbsp;
* S-box

  is a basic component of symmetric key algorithms which performs substitution. In block ciphers,\
  they are typically used to obscure the relationship between the key and the ciphertext — Shannon's\
  property of confusion. In general, an S-box takes some number of input bits, m, and transforms them\
  into some number of output bits, n, where n is not necessarily equal to m. An m×n S-box can be\
  implemented as a lookup table with 2m words of n bits each. Fixed tables are normally used, as in\
  the Data Encryption Standard (DES), but in some ciphers the tables are generated dynamically from\
  the key (e.g. the Blowfish and the Twofish encryption algorithms).\
&nbsp;
&nbsp;
* meet in middle attack

  is a generic space–time tradeoff cryptographic attack against encryption schemes that rely on\
  performing multiple encryption operations in sequence.\
&nbsp;
&nbsp;
* cryptography

  is the practice and study of techniques for secure communication in the presence of third parties\
  called adversaries. More generally, cryptography is about constructing and analyzing protocols that\
  prevent third parties or the public from reading private messages; various aspects in information\
  security such as data confidentiality, data integrity, authentication, and non-repudiation are\
  central to modern cryptography. Modern cryptography exists at the intersection of the disciplines of\
  mathematics, computer science, electrical engineering, communication science, and physics.\
  Applications of cryptography include electronic commerce, chip-based payment cards, digital\
  currencies, computer passwords, and military communications.\
&nbsp;
&nbsp;
* cryptology

  is a branch of mathematics which deals with both cryptography and cryptoanalysis.\
&nbsp;
&nbsp;
* cryptographic protocol

  is an abstract or concrete protocol that performs a security-related function and applies\
  cryptographic methods, often as sequences of cryptographic primitives. A protocol describes how the\
  algorithms should be used. A sufficiently detailed protocol includes details about data structures\
  and representations, at which point it can be used to implement multiple, interoperable versions of\
  a program.
&nbsp;
&nbsp;
* encryption

  is the process of encoding information. This process converts the original representation of the\
  information, known as plaintext, into an alternative form known as ciphertext. Ideally, only\
  authorized parties can decipher a ciphertext back to plaintext and access the original information.\
  Encryption does not itself prevent interference but denies the intelligible content to a would-be\
  interceptor. For technical reasons, an encryption scheme usually uses a pseudo-random encryption key\
  generated by an algorithm. It is possible to decrypt the message without possessing the key, but, for\
  a well-designed encryption scheme, considerable computational resources and skills are required. An\
  authorized recipient can easily decrypt the message with the key provided by the originator to\
  recipients but not to unauthorized users. Historically, various forms of encryption have been used\
  to aid in cryptography. Early encryption techniques were often utilized in military messaging. Since\
  then, new techniques have emerged and become commonplace in all areas of modern computing. Modern\
  encryption schemes utilize the concepts of public-key and symmetric-key. Modern encryption techniques\
  ensure security because modern computers are inefficient at cracking the encryption.\
&nbsp;
&nbsp;
* decryption

  is the reverse, in other words, moving from the unintelligible ciphertext back to plaintext.\
&nbsp;
&nbsp;
* DES (Data Encryption Standard)

  is a symmetric-key algorithm for the encryption of digital data. Although its short key length of\
  56 bits makes it too insecure for applications, it has been highly influential in the advancement\
  of cryptography. DES has been remarkably resistant to cryptanalysis, but its short key length makes\
  it vulnerable to a brute-force attack where all possible keys are tried one by one until the correct\
  key in found.\
&nbsp;
&nbsp;
* RSA (Rivest-Shamir-Adleman)

  is a public-key cryptosystem that is widely used for secure data transmission. It is also one of the\
  oldest. In a public-key cryptosystem, the encryption key is public and distinct from the decryption\
  key, which is kept secret (private). An RSA user creates and publishes a public key based on two\
  large prime numbers, along with an auxiliary value. The prime numbers are kept secret. Messages\
  can be encrypted by anyone, via the public key, but can only be decoded by someone who knows\
  the prime numbers.\
&nbsp;
&nbsp;
* digital signature

  is a mathematical scheme for verifying the authenticity of digital messages or documents. A valid\
  digital signature, where the prerequisites are satisfied, gives a recipient very strong reason\
  to believe that the message was created by a known sender (authentication), and that the message\
  was not altered in transit (integrity). Digital signatures employ asymmetric cryptography. In many\
  instances they provide a layer of validation and security to messages sent through a non-secure\
  channel: Properly implemented, a digital signature gives the receiver reason to believe the message\
  was sent by the claimed sender\
&nbsp;
&nbsp;
* steganography

  is the practice of concealing a file, message, image, or video within another file, message, image,\
  or video.\
&nbsp;
&nbsp;
* cryptographic algorithm

  is a mathematical function which uses plaintext as the input and produces ciphertext as the output\
  and vice versa.\
&nbsp;
&nbsp;
* **distinguisher evaluation in BoolTest**

  An empirical test of randomness consists (in general) of the following steps:\
  1. Compute the histogram H of some features (within data).
  2. Compute (transform the histogram to) the observed test statistic s_obs.
  3. Compute the null distribution (exact or its close approximation) D(x)of the test statistic under\
     the null hypothesis (random data).
  4. Compute the p-value from s_obs using the distribution D(x).\

  In our approach, the histogram of results of the boolean function f(x1, ···, xm) of m variables\
  applied to non-overlapping m-bit blocks of the sequence is computed. Our test statistic is Z-score\
  (Sheskin, 2003) defined as: Z-score = (#1−pn) / (√p(1−p)n), which normalize a binomial distribution\
  B(n,p). Binomially distributed variable #1 is normalized to Z-score which is distributed normally.\
  P-value can be directly computed from the Z-score. Figure 1b illustrates the relation of\
  a Z-score (standardly expressed in the units of standard deviation x·σ) and the corresponding\
  p-value (area of two tails).\
  The symbol p denotes the probability that the result of boolean function f is equalto 1 for random\
  input. The symbol n denotes the number of non-overlapping blocks (of m bits) in the analyzed sequence\
  (not the number of bits). Similarly, as in the Monobit test, our histogram consists of two frequencies\
  #1 and #0, but only #1 is computed (#0 = n − #1) and is used for the evaluation. The only difference\
  is that the expected probability p is not p = 0.5. In general, p is arbitrary value from the interval\
  [0, 1] which depends on the given boolean function f.\
&nbsp;
&nbsp;
