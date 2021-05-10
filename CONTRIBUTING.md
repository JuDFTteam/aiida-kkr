# Contributing guide of aiida-kkr

## Introduction

First off all, thank you for considering contributing to aiida-kkr.

We are always looking for helping hands. There are many ways you can contribute from writing tutorials and improving the documentation over bug fixes to the development of new features (e.g. workflows).

## Ground Rules

* Create issues for any major changes and enhancements that you wish to make. Discuss things transparently and get community feedback.
* Keep feature versions as small as possible, preferably one new feature per version.
* Write tests for new features you implement.
* Document the source code.
* Be welcoming to newcomers and encourage diverse new contributors from all backgrounds. See the [Python Community Code of Conduct](https://www.python.org/psf/codeofconduct/).

## Branch structure of aiida-kkr

We try to follow the gitflow model of branching where the latest stable version is on the `master` branch.
This is the branch from where releases are created.

In addition we have a `develop` branch where current developments are collected / done. For developing a new feature that might take some time to implement please consider creating a new branch.

## How to report a bug

Any security issues should be submitted directly to p.ruessmann@fz-juelich.de.
In order to determine whether you are dealing with a security issue, ask yourself these two questions:
* Can I access something that's not mine, or something I shouldn't have access to?
* Can I disable something for other people?

If the answer to either of those two questions are "yes", then you're probably dealing with a security issue. Note that even if you answer "no" to both questions, you may still be dealing with a security issue, so if you're unsure, just email us at p.ruessmann@fz-juelich.de.


When filing an issue, make sure to answer these five questions:

1. What version of python are you using?
2. What operating system and processor architecture are you using?
3. What did you do?
4. What did you expect to see?
5. What did you see instead?
