# Git workflow guidelines

- Never work on the `master` branch. Create a feature branch, work there, then merge into `master`.
- Try to keep feature braches short-lived, or at least keep them up to date with `master`. This can save you a lot of headaches when merging into `master`.
- When merging feature branches:
    1. **Announce** what you're doing on the Slack channel. This is to avoid race conditions where multiple people are attempting to rebase and merge at the same time.
    2. **Pull** changes to the master branch from the central repository.
    3. **Clean up** your feature branch and **update it** with changes from `master`. Two ways of doing this:
        - Rebasing `feature` onto `master`: `git checkout feature && git rebase -i master` (give [this article](https://www.atlassian.com/git/tutorials/merging-vs-rebasing) a read if you're unfamiliar with the **golden rule of rebasing**!).
        - Merging `master` into `feature`: `git checkout feature && git merge master`. This has the upside of being safer if you're unsure about rebasing, but creates a messier history.
    4. **Check** the resulting codebase to make sure nothing is broken.
    5. **Merge** `feature` into `master`; at this point it should be a simple fast-forward commit with no conflicts to resolve: `git checkout master && git merge feature`.
