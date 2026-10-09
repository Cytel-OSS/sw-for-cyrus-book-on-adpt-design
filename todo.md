## To check with Cyrus
- "Chapter 5/Sec 5.2.2/MinDeltaDesignSimulations-PD.r" - does not produce any plots, but performs some kind of simulation. Do we need to upload it to GitHub?


## Changes to be made to R code
At the top of the R script, between `rm(list = ls())` and `library()` statements, add the following:
```
# To prevent R from creating an Rplots.pdf file when this script is executed non-interactively.
options(device = function(...) {
  grDevices::pdf(NULL)
})
```

Then, after the `library()` statements, add the following code to source `helper.R` file:
```
## Paths and directories
# Sourcing the helper functions from the tools directory
tryCatch(
  source(file.path("tools", "helper.R"), local = .GlobalEnv),
  error = function(first.error) {
    source(
      file.path("..", "..", "tools", "helper.R"),
      local = .GlobalEnv
    )
  }
)
```

After that, call the function `get.base.dir()` from `helper.R` to set `base.dir` as shown below:
IMPORTANT: Modify `sec.dir` var appropriately based on your section folder's relative path.
```
# Setting up the base directory for executing this code
base.dir <- get.base.dir(
  cur.dir = getwd(),
  sec.dir = "Chapter 5/Sec 5.2.1/"
)
```

Change Apurva's existing code for setting base.dir. She was working with a different folder structure with different input/output folders, which was unnecessarily complex and was not portable.
I am changing this to assume that the R code and any datafile it requires reside in the same folder (the section folder), and the R code will create the output in that same folder.

Finally, at the bottom of the code, before the `ggsave()` call, remove Apurva's existing assignment of `fig.dir` and set it to `base.dir` so that the output figure also gets created in the same section folder:
```
fig.dir <- base.dir
```

## Notes from 22 Sep 2026 meeting
I reviewed the full **22 September 2026 meeting transcript**. Unlike the previous two meetings, this session was **not primarily about AdaptGMCP implementation**. Cyrus explicitly said the purpose was to organize the software/resources accompanying his upcoming book, particularly the examples in **Chapters 4 and 5**, and to make those resources usable by readers. AdaptGMCP Meeting - 22Sep2026

### Executive summary

The central objective was to turn Cyrus's existing collection of R programs, East workbooks, and generated results into a **portable, organized, documented software resource that can be distributed through GitHub** and linked clearly to the relevant sections of the book.

The immediate problem is that the current material was developed on Cyrus's machine and contains **hard-coded paths and assumptions about folder structure**, so simply copying the individual files does not work. The agreed approach is to preserve the required folder structures, distribute them as ZIP files initially, and then reorganize them in GitHub using a structure aligned with **book sections rather than arbitrary filenames**. AdaptGMCP Meeting - 22Sep2026

The meeting also established a sensible documentation philosophy: readers will read the relevant chapter, go to the corresponding software folder, follow a short "how to run this" guide, and obtain the figures/tables discussed in the chapter. The software should be sufficiently self-contained that a reader can run an example without having to reverse-engineer the author's development environment. AdaptGMCP Meeting - 22Sep2026 AdaptGMCP Meeting - 22Sep2026

The **Mehta-Pocock example** was identified as the easiest place to start. More complicated JT, pills, and Bayesian examples will be tackled progressively later. AdaptGMCP Meeting - 22Sep2026

### Important points

**1. This meeting was about book software resources, not AdaptGMCP**

Cyrus explained that the immediate need is to make the book usable for readers who need access to statistical software. The required software environment includes **R and East**, with East already having established documentation. The focus of this work is therefore primarily on making the **R-based examples accessible and portable**. AdaptGMCP Meeting - 22Sep2026

**2. Existing code is organized into sets corresponding to Chapter 5**

The examples currently exist in several sets. Ani reviewed them against the chapter and identified that some later examples, particularly around Section 5.5, were difficult to execute because they depend on a particular directory structure and hard-coded paths. AdaptGMCP Meeting - 22Sep2026

Cyrus explained that the intended design was actually to keep **code and results in separate folders**, so that the whole structure can be moved together and executed elsewhere. AdaptGMCP Meeting - 22Sep2026

**3. Portability is a key requirement**

Ani proposed modifying the programs so that the required folder structure can be reproduced on another machine and the whole resource can eventually be put on GitHub. The goal is that someone downloading the repository should not have to discover Cyrus's original directory structure manually. AdaptGMCP Meeting - 22Sep2026

Cyrus agreed and emphasized that the examples should actually be tested after this conversion so that they run and produce the intended outputs. AdaptGMCP Meeting - 22Sep2026

**4. East examples and R examples should be treated differently**

A useful distinction emerged:

- **East-based examples:** East already has good documentation, so the main requirement is helping readers obtain access to East.
- **R-based examples:** These need the accompanying code, organization and instructions because the software itself is not being distributed as a turnkey GUI application. AdaptGMCP Meeting - 22Sep2026

**5. The book and software should be explicitly linked**

Cyrus wants each software collection to correspond directly to a section of Chapter 5. The idea is that a reader who reaches a particular section can immediately find the software associated with that section and use it to reproduce or explore the material. AdaptGMCP Meeting - 22Sep2026

The proposed organization is therefore something like:

> Chapter 5 → Section 5.x → corresponding software folder → code/results/instructions.

Rather than using long descriptive folder names everywhere, Ani suggested short section-based folder names with a parent directory and an **index/README** explaining what each folder contains. Cyrus agreed that this could work well. AdaptGMCP Meeting - 22Sep2026

**6. A concise "how to run" document is important**

For each example, there should be a short document—perhaps half a page or one page—that tells the reader:

- what the example does,
- which program to run first,
- what subsequent programs to run,
- what output to expect, and
- what prerequisites/packages are needed.

Ani specifically suggested documenting required R packages so that users do not discover dependencies only through runtime errors. AdaptGMCP Meeting - 22Sep2026

**7. Some examples involve multiple sequential R programs**

The more complex examples are not necessarily one-script examples. Some consist of several R programs where one program produces intermediate output, often a CSV file, that is then consumed by another program. The user documentation therefore needs to explain the execution sequence clearly. AdaptGMCP Meeting - 22Sep2026 AdaptGMCP Meeting - 22Sep2026

**8. Some computations are intentionally intensive**

One of the JT/pills design examples first performs an intensive calculation and saves its results to CSV files; subsequent code uses those saved results. Cyrus explicitly wants users warned that the initial calculation can take a long time. AdaptGMCP Meeting - 22Sep2026

For users interested only in the JT designs, the code can be stopped earlier, rather than running the paired-design calculations as well. AdaptGMCP Meeting - 22Sep2026

**9. The examples should minimize interactive steps where possible**

For the more complicated sets, Cyrus and Ani discussed determining whether the programs can be made to run sequentially rather than requiring interactive intervention. The preference is to automate as much as possible so that a reader can follow the documented sequence and obtain the required figures. AdaptGMCP Meeting - 22Sep2026

**10. Gradual implementation was explicitly preferred**

The software resource collection does not need to be completed all at once. Cyrus suggested doing the examples incrementally over the coming months, starting with the easier ones and progressively addressing the more complicated examples. AdaptGMCP Meeting - 22Sep2026

### Agreed action items

| Owner | Action |
|---|---|
| **Ani** | Review the supplied R programs and **remove hard-coded machine-specific paths**, making the examples portable. AdaptGMCP Meeting - 22Sep2026 |
| **Ani** | Organize the software into a **GitHub-friendly folder structure aligned with Chapter 5 sections** rather than the current arbitrary structure. AdaptGMCP Meeting - 22Sep2026 |
| **Cyrus** | Provide the required source material as **ZIP files preserving the original directory structures**, rather than sending individual files. AdaptGMCP Meeting - 22Sep2026 |
| **Cyrus** | Send the required ZIP(s) through **Teams** because the files are too large / ZIP attachments may be blocked by email security controls. AdaptGMCP Meeting - 22Sep2026 |
| **Ani** | Start by working on the **Mehta-Pocock example**, identified as the easiest initial example. Then proceed to the next example(s). AdaptGMCP Meeting - 22Sep2026 |
| **Ani** | Ensure each portable example is actually **executed and verified** after modification. AdaptGMCP Meeting - 22Sep2026 |
| **Ani** | Add concise documentation/README material describing **what each folder contains, how to run it, expected output, and required R packages**. AdaptGMCP Meeting - 22Sep2026 |
| **Cyrus** | Gradually add/adjust explanatory comments and descriptions so the examples are easier for readers to understand and correspond clearly to the book. AdaptGMCP Meeting - 22Sep2026 |
| **Cyrus + Ani** | Work through the more complex multi-step examples progressively, determining which can be made non-interactive and documenting the required execution sequence. AdaptGMCP Meeting - 22Sep2026 |

### Resulting direction

The overall design that emerged is quite clear:

**Book section → GitHub software folder → short README/instructions → R programs + supporting files/results → reproducible figure/table.**

This should make the software a genuine companion to the book rather than simply a collection of source files. Cyrus's closing observation was that the chapters are technically difficult, and having runnable software alongside them should make the material considerably easier for readers to understand. AdaptGMCP Meeting - 22Sep2026

The next AdaptGMCP meeting was expected to return to the package work the following week. AdaptGMCP Meeting - 22Sep2026
