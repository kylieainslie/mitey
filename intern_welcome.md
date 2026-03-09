# Welcome to the Internship
**Project:** Extending the {mitey} R Package
**Duration:** 4 months (16 weeks)

---

## What This Internship Is About

This internship has two goals that should feel connected throughout: learning R in the context of a software project, and performing data analysis on infectious disease data. The project you're contributing to — the [`{mitey}` R package](https://kylieainslie.github.io/mitey/) — is used to estimate how quickly infectious diseases spread between people. Your contribution will make it more flexible and useful for analysing a wider range of outbreaks.

By the end of the internship you'll have written an extension to `{mitey}` that extends a published statistical method, applied it to real outbreak data, and produced a vignette and short write-up documenting what you found.

---

## Your Role

**What you will build:**

`si_estim_flex()` — a new function that extends the existing `{mitey}` estimation method to support a flexible number of mixture components (rather than a fixed four). This is the core technical contribution of the internship, and the Methodological Development Plan is your technical guide for building it.

**What you'll do with the tools:**

Apply `si_estim_flex()` to a real historical outbreak dataset. Fit models with different numbers of components, compare the fits, interpret the results in epidemiological terms, and write it up as a package vignette.

**Final deliverables:**
- A working `si_estim_flex()` function for K=2–4 components
- One real data analysis vignette included in the package
- A short technical write-up: background, methods, results, what you learned
- A presentation to the team

---

## Four-Month Overview

| Month | Theme | Goal by end of month |
|-------|-------|----------------------|
| 1 | Onboarding, R learning, codebase understanding | Comfortable writing R; understand how {mitey} works; prototype scaffold started |
| 2 | Core implementation | `si_estim_flex()` working and tested for K=2–4; passing R CMD check |
| 3 | Validation + data exploration | Implementation validated; real dataset explored; ready to write up |
| 4 | Data analysis + write-up | Vignette written; report and presentation delivered |

The plan is designed to flex based on how things go. Monthly check-ins with your supervisor are built in to adjust the focus as needed — so don't worry if things take longer than expected in the early weeks. The R learning phase is deliberately unhurried.

---

## Background Reading

Work through these alongside the technical setup in Weeks 1–2. They're listed in the order that makes most sense to read them.

### Epidemiology — start here

Before looking at any code, it's worth understanding what the package is doing and why. The key concepts are:

**Serial interval** — the time between symptom onset in a person who infects someone else, and symptom onset in the person they infected. It's a measure of how fast a disease spreads and it differs from the *incubation period* (time from infection to symptoms in one person) and the *generation time* (time between infections, which can't be directly observed).

**ICC intervals** — "index case-to-case" intervals are the raw data {mitey} uses. In a household outbreak, you observe the times between cases — but some of those cases infected each other directly, while others were both infected by a common source. The mixture model is the tool for separating those signals.

1. [**Vink et al. (2014)**, "Serial intervals of respiratory infectious diseases: a systematic review and analysis"](https://doi.org/10.1093/aje/kwu209) (*American Journal of Epidemiology*) — the paper that {mitey} implements. Read this first. Focus on the biological motivation and the intuition for why a mixture model is needed; the methods section will make more sense after you've worked through the statistical reading below.

2. [**Svensson (2007)**, "A note on generation times in epidemic models"](https://doi.org/10.1016/j.mbs.2006.10.010) (*Mathematical Biosciences*) — short and clear. Explains exactly why serial interval, generation time, and incubation period are different things — and why serial intervals can even be negative.

3. [**Nishiura (2007)**, "Time variations in the transmissibility of pandemic influenza in Prussia, Germany, from 1918–19"](https://doi.org/10.1186/1742-4682-4-20) (*Theoretical Biology and Medical Modelling*, open access) — a worked example of the kind of analysis you'll be doing. Useful for seeing the concepts applied to a real outbreak.

4. [**Lessler et al. (2009)**, "Incubation periods of acute respiratory viral infections: a systematic review"](https://doi.org/10.1016/S1473-3099(09)70069-6) (*Lancet Infectious Diseases*) — useful background on why parametric distributions (Normal, Gamma, lognormal) are commonly used to model infectious disease timing data.

5. [**Anderson & May (1991)**, *Infectious Diseases of Humans: Dynamics and Control*](https://global.oup.com/academic/product/infectious-diseases-of-humans-9780198540403) (Oxford University Press) — not required, but your supervisor may point you to specific chapters for broader context.

### Statistics — read in Weeks 2–3

These build the statistical foundation you'll need before implementing the EM algorithm. The sequencing matters: read the epidemiology papers first so the statistics has a concrete problem to attach to.

**Maximum likelihood estimation — start here**

The concept of a likelihood function and what it means to maximise it is the single most important statistical idea for this project. Make sure you understand this before looking at any {mitey} source code.

- [**StatQuest: Maximum Likelihood, clearly explained**](https://www.youtube.com/watch?v=XepXtl9YKwc) (YouTube, ~10 min) — watch this first. Visual and concrete, no prerequisites.
- [**StatQuest: Probability vs Likelihood**](https://www.youtube.com/watch?v=pYxNSUDSFH4) (YouTube, ~5 min) — clears up a distinction that trips up almost everyone coming from an engineering background.

**The EM algorithm**

- [**Do & Batzoglou (2008)**, "What is the expectation maximization algorithm?"](https://doi.org/10.1038/nbt1406) (*Nature Biotechnology*) — three pages, written explicitly for practitioners without a deep statistics background. Works through a concrete coin-flipping example. Read this before looking at the {mitey} source code.
- [**StatQuest: The EM Algorithm, clearly explained**](https://www.youtube.com/watch?v=REypj2sy_5U) (YouTube, ~20 min) — good visual companion to Do & Batzoglou.

**Mixture models**

- [**Reynolds (2009)**, "Gaussian Mixture Models"](https://doi.org/10.1007/978-0-387-73003-5_196) (*Encyclopedia of Biometrics*) — concise overview of Gaussian mixture models and the EM algorithm as applied to them. Read after Do & Batzoglou.
- **Bishop (2006)**, *Pattern Recognition and Machine Learning*, Chapter 9 — "Mixture Models and EM". More detailed treatment for when you want depth during implementation. The [full PDF is freely available](https://www.microsoft.com/en-us/research/uploads/prod/2006/01/Bishop-Pattern-Recognition-and-Machine-Learning-2006.pdf).

**Suggested reading order:**
Vink et al. → Svensson → StatQuest MLE videos → Do & Batzoglou → Reynolds → {mitey} source code

---

### Month 1 — Onboarding, R Learning, Codebase Understanding

**End-of-month milestone:** R proficiency for package development; understands the existing {mitey} EM algorithm; has begun prototype scaffold for `si_estim_flex()`

---

#### Week 1 — Institutional Onboarding & R Fundamentals

**Institutional onboarding (first 1–2 days — complete before anything else)**
- [ ] Collect ID card and building access
- [ ] Complete occupational health and safety induction (mandatory before starting lab or office work)
- [ ] Complete any required institutional trainings (data privacy, research ethics, IT security, or similar — confirm list with supervisor on Day 1)
- [ ] Set up institutional accounts: email, VPN, shared drives, any required systems
- [ ] Review internship goals, schedule, communication norms, and weekly meeting time with supervisor

**Technical setup**
- [ ] Install and configure: 
  - [ ] [R](https://cran.r-project.org/)
  - [ ] [RStudio](https://posit.co/downloads/) or [Positron](https://positron.posit.co/download.html)
  - [ ] [Git](https://git-scm.com/install/); 
  - [ ] connect Git to [GitHub](https://github.com/) (setup a GitHub account if you do not already have one) 
  - [ ] clone [`{mitey}` repo](https://github.com/kylieainslie/mitey) and check out your working branch:
    ```bash
    git clone https://github.com/kylieainslie/mitey.git
    cd mitey
    git checkout simon
    ```


**R and project orientation**
- [ ] R basics: vectors, lists, data frames, indexing — compare to Python lists/dicts/pandas
- [ ] R vectorized operations (vs. NumPy/pandas broadcasting); writing short scripts with functions, loops, conditionals
- [ ] Explore R's probability distribution functions: `dnorm`, `dgamma`, `pnorm`, `rnorm`, etc.
- [ ] Install and run [Swirl](https://swirlstats.com/students.html) ("R Programming" course) for interactive syntax practice
- [ ] Read [{mitey} README](https://kylieainslie.github.io/mitey/) and run existing examples; begin personal notes document

---
