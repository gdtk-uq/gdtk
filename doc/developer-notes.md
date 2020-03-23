# Notes on the development process.

Contents:

* Introduction
* Running the automated smoke tests
* Making a commit
* What to commit and what not to commit
* Odds and ends

## Introduction

At the moment, there is essentially a single branch in the main repository. 
This is because I like a linear history, or something close to it, 
for code that is relatively young, has no long-lived versions and 
is worked on by a fairly small group of people. 
As the project grows, I'm sure that we'll end up with long-term branches, 
at least for code revisions that people expect to be able to cite, 
say Eilmer4 version 4.0.0. 
Whenever there are multiple branches, someone needs to maintain them 
and merge code from one to the other. 
We don't wish to pay that price, yet.

So, given that we have a linear history of the master branch, 
what's the best way to develop against it?

I like to work on a single issue in any one development session. 
Starting with a code base that is in a known good state, 
I like to check for incoming commits (from other people) and 
update to the latest, if there are some incoming commits. 
Look at the new changes and see if it affects me or my work 
for the present session. 
If all is good, go ahead and make code changes 
related to a particular issue. 
Sometimes the issue involves a lot of new code or edits. 
Sometimes, it might be just fixing an irritating bit of trivia; 
it's not good to leave stuff broken even if it is small and 
seemingly trivial, like badly worded comments with speling errors.

## Running the automated smoke tests

When I'm done, I confirm that the code is in a working state 
for my new stuff and then run the automated test set. 
Typically I use the commands:

    $ cd ~/work/eilmer4-test/
    $ rsync -av --delete ~/dgd/examples/eilmer/ .
    $ ./eilmer-test.tcl

Have dinner, or at least a cup of tea, while waiting about an hour 
for the full set to complete on my.

Fix any tests that complain and retest.

Once, I'm happy that things are working well and 
that I haven't broken any other code, I'm ready to commit.

## Version Control

Eilmer uses the venerable distributed version control system git to manage its repository of source code.
This repository is a hash-tree of changes hidden in the \.git directory next to the source files, which tracks every change made to every line of code, going back to the beginning of the project.
Version control is critical to managing backups, finding bugs, and integrating work from multiple developers, so developing the code requires a bit of knowledge about how it works in addition to just knowing the terminal commands.
If you're unfamiliar with git, it is a good idea to check out a tutorial such as https://www.atlassian.com/git, to get familiar with the basic concepts.
The next few sections document outlines the basic recommended workflow for making changes to Eilmer.

## Developing with git

Before beginning any work, it is a good idea to start with a fresh pull of the repository.
From inside an old dgd directory type:

    $ git pull

pull is a multipart command that checks the server where Eilmer is hosted online for new changes, integrates those changes into your local repo, and changes the source code in your dgd directory to match the newest version of each file.

At this point you're free to make changes to the code, add new features, debug old issues, etc.
git has many tools that might help you along the way, such as:

    $ git status # shows files changed since last commit
    $ git diff # shows differences between your current files and last commit
    $ git log --graph --all # shows recent commits and history graph
    $ git diff VERSIONHASH1 VERSIONHASH2 # shows changes between old commits
    $ git stash # saves uncommitted changes so you can go back to the last commit
    $ git stash pop # loads uncommitted changes from git stash
    $ git checkout -- FILENAME # reverts FILENAME to last commit

Once you've finished your changes and tested that they work, it is a good idea to run the automated tests in examples/eilmer:

    $ cd examples/
    $ cp -r eilmer eilmer.test
    $ cd eilmer.test
    $ ./eilmer-test.rb

Running these tests requires all of the versions of eilmer (including the steady-state and complex number versions), as well as some extra dependencies such as python-sympy and ruby, and should take about a hour-and-a-half on a good machine.
Once you're satisfied that your changes have not harmed the code in some way, it's time to commit them to your local copy of the repository, and push it back out to the server so other people can use your changes.

## Adding Changes to the repository

To start with, it is a good idea to run make clean and remove untracked files.
This will clean up the output of git status and help you keep track of what's going on.

    $ cd src/eilmer
    $ make clean
    $ git status

The main issue with pushing changes to the repository is that someone else may have added commits there since you started working.
To check for incoming changes use a git pull:

    $ git pull

If successful, this will add new changes into your repository and update the source code files to match them.
If this is the case then skip ahead to the "Making a commit" section.
However, if you have made changes to a file that is also set to be changed, git pull will abort and explain politely that it cannot continue.
If this is the case you will need to perform a merge.

## Merging

Modern version control systems contain clever routines for automatically merging the changes in a file that has a diverging development history.
The easiest way to perform a merge in git without making a commit is to stash your changes, pull the new changes and unstash them:

    $ git stash
    $ git pull
    $ git stash pop

Under the hood, git stash takes local changes you've made to the source code and stashes them into a temporary commit that is not added to the main development branch.
If you were to check any of the files at that point you would find that they have been reset to how they were before you started working on them.
Because of this, git pull will now complete, updating the files to bring in new changes.
At this point git stash pop will automatically merge the conflicting changes together, slicing in the lines of code you changed into each file, using their shared ancestor as a reference point.
 
Most commonly this merge will succeed automatically.
However if your changes were to the same part of the same files you will get a merge conflict and the automerge will stop.
Solving a merge conflict requires manually deciding which line of code should win at each contested point, by opening up the file and going to the contested lines and editing them yourself.

Once the merge has been completed it may be a good idea to run the automated tests again, to make sure the automatic system has not broken the code.
This is especially true in the case of a conflicted or messy merge, or if your changes were to an important part of the code that interacts with many other parts.
Keeping up to date with the main repository and doing small frequent commits should reduce the risk of these issues, but sometimes they are unavoidable.

## Making a commit

When merging is complete it is time to commit the changes to your local copy of the repository.
Take a look at git status to make sure that the changes are what you expect and then add the files to staging area and commit them.
If you have added a new file to the code make sure to add it too, or it won't show up in other people's copies of the code!

    $ git add file1.d file2.d newfile.d
    $ git commit -m "A commit message"

Make sure to add a good descriptive commit message, this will help other people understand what you have done, and help you in case you have to come back to this commit and debug anything.

At this point your repository should be in a state where it can be pushed to the master copy other people will pull down and use.
A simple command to do this is:

    $ git push

The server will ask you for credentials, and if you are authorised, will attempt to push your changes.
Rarely, someone may have performed a push since your last pull command, perhaps while you were doing a second round of tests.
In this case git will prevent you from pushing and tell you to pull again.
If this happens you will have to pull again and fix a divergent commit history, using either git merge or git rebase, then try git push.
The best way to avoid this situation is to minimise the amount of time between git pull and git push.
In the case of a complicated or messy merge that takes some time, simply run git pull again before committing to ensure that everything is up to date.

## What to commit and what not to commit

### Yes

* A single feature or fix, possibly over several files.
At the other end of the scale, maybe only a character or two.
* Manually entered scripts and data sets.  
Files that you have laboured over and don't want to lose.

### No

* Complicated change sets across many files and with multiple, 
unrelated features or fixes.  It is usually desirable to split 
these into multiple commits with one issue per commit.
* Source code changes that are just white-space changes made because
your editor behaves differently to `emacs` on Linux 
with `c-basic-offset` set to 4 and `indent-tabs-mode` set to `nil`.
Note that Lua files have 3 spaces per indent level.
Gratuitous changes to whitespace add to the burden of reviewing commits.
* Program executable codes that can be regenerated from the source.
Do a `make clean` before committing.
* Other binary files such as gzipped solution files or PDF files,
unless they are "golden solutions" that will never change.
Binary files rapidly add bloat to the repository.
* Backup files produced by your editor.  Use a `.hgignore` file in your
working tree.

## Odds and ends

### Emacs
Here is the content of my `.emacs` file:

    peterj@helmholtz ~ $ cat .emacs
    (custom-set-variables
      ;; custom-set-variables was added by Custom.
      ;; If you edit it by hand, you could mess it up, so be careful.
      ;; Your init file should contain only one such instance.
      ;; If there is more than one, they won't work right.
     '(c-basic-offset 4)
     '(indent-tabs-mode nil))
    (custom-set-faces
      ;; custom-set-faces was added by Custom.
      ;; If you edit it by hand, you could mess it up, so be careful.
      ;; Your init file should contain only one such instance.
      ;; If there is more than one, they won't work right.
     )
    (autoload 'd-mode "d-mode" "Major mode for editing D code." t)
    (add-to-list 'auto-mode-alist '("\\.d[i]?\\'" . d-mode))
    
    (autoload 'lua-mode "lua-mode" "Lua editing mode." t)
    (add-to-list 'auto-mode-alist '("\\.lua$" . lua-mode))
    (add-to-list 'interpreter-mode-alist '("lua" . lua-mode))

### Git configuration
I did the following:

    peterj@helmholtz ~ $ git config --global user.name "Peter Jacobs"
    peterj@helmholtz ~ $ git config --global user.email "peterj@mech.uq.edu.au"

Note that the user.name is very explicit: it leaves no doubt as
to who to blame!
On a more serious note, please use your full name for user.name.
This will help identifying the contributors to the code over 
the long-term history of development.

### Other editors
If you are using an editor other than emacs, please configure your editor 
so that it uses an indentation level of 4 spaces for D-language files
and 3 spaces for Lua files.
Use spaces for indentation, not tab characters.
Although we have moved away from using tab characters, 
we continue with the convention that tab stops are every 8th column.

An easy way to configure your editor of choice
is to use the plugins provided at:

http://editorconfig.org/

In the same directory as this note is an `.editorconfig` that
has the appropriate settings to play nicely with our preferred
indentation style. The file is named `editorconfig` here so
that it isn't hidden from view, but it should be renamed as
`.editorconfig` if you want to use it with a plugin.
