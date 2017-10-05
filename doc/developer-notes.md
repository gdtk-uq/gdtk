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

So, given that we have a linear history, 
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

Have dinner, or at least a cup of tea, while waiting about half an hour 
for the full set to complete on my slow workstation. 
On a modern system the tests might take only 15 minutes but 
then I wouldn't get a good chance to relax with dinner.

Fix any tests that complain and retest.

Once, I'm happy that things are working well and 
that I haven't broken any other code, I'm ready to commit.

## Making a commit

Before making the commit, I look to see if there are any other commits 
that have appeared since the start of my session and, 
if there are any, pull them.
I usually also update but there is a finite probability that 
the other commits will interfere with my new code.
If that happens, the merge tool usually does a good job 
of integrating the pulled changes into my current tree of files.
Often, there are no files to be merged and the update proceeds cleanly 
but I now need to check that code is still working.
Sometimes I get a second cup of tea out of this, 
as the automated tests run again.
Occasionally, I decide to commit my new stuff immediately and 
then run the tests to check that my assumption that all is good 
is really true.

So, I make my commit on top of the full history 
of the bitbucket repository, thus keeping a linear history.
If you browse the actual revision history, you will see 
various small branches followed by merges where 
different developers have made commits in parallel.
Sometimes, if the parallel work has been done on distinct files, 
the merge is trivially done by the Mercurial revision control system.
Sometimes edits collide and someone has to manually check/do the merge.
It's nice to avoid this extra work.

Once you have your commit with your new work sitting at the tip of your 
copy of the repository, it's time to push it back to 
the bitbucket hosted repository.
There is a race condition at this point because other people will 
proceed to do their development work at their own pace and 
may push back their commits before you push yours.
This occasionally happens but usually is not a big deal.
Pull and update, and possibly merge again, then try to push back.
So far, I've not needed to do this merge cycle more than twice, 
but we'll see how things go as more people participate in the code development.

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
your editor behaves differently to `emacs` with `c-basic-offset` set to 4.
Note that Lua files have 3 spaces per indent level.
Tabs are fine so long as they are set to a width of 8 spaces.
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
     '(c-basic-offset 4))
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

### Mercurial

Here is the content of my `.hgrc` file:

    peterj@helmholtz ~ $ cat .hgrc 
    [ui]
    username = Peter Jacobs <peterj@mech.uq.edu.au>
    
    [extensions]
    hgext.convert =

Note that the 'username' is very explicit: it leaves no doubt as
to who to blame!
On a more serious note, please use your full name in the 'username' field.
This will help identifying the contributors to the code over 
the long-term history of development.

### Other editors
If you are using an editor other than emacs,
please configure your editor so that it uses 8-spaces
per tab-stop and indentation level of 4.
An easy way to configure your editor of choice
is to use the plugins provided at:

http://editorconfig.org/

In the same directory as this note is an `.editorconfig` that
has the appropriate settings to play nicely with our preffered
indentation style. The file is named `editorconfig` here so
that it isn't hidden from view, but it should be renamed as
`.editorconfig` if one wants to use it with a plugin.