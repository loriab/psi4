Let me try to document the procedure (besides the actual version bump, which is already documented). This issue can be scratch space, so everyone can post easily.

Update copyright year in licenses:
user@host:~/path/to/psi4 $ grep -rl "(c) 2007-2017" * | xargs sed -i '' "s/(c) 2007-2017/(c) 2007-2018/g"

on copyright, psi4/ and docs/ dirs need to be done.
on Linux, drop the '' in above command

* run `make sphinxman` at least once by hand and check in resulting psifiles.py and all the updated and new `samples/` files and dirs.

* figure out any new `Additional Contributors` since last release and edit header.py accordingly. Get their particulars for codemeta.json. invite any >=3 PRs to GH org.

* follow the procedure outlined in [cbcy](https://github.com/psi4/psi4meta/blob/master/conda-recipes/conda_build_config.yaml) and [poodle](https://github.com/psi4/psi4meta/blob/master/psinet-nightly/kitandkapoodle.py) (I haven't pushed the versions with enhanced comments) to
  * update upstream and downstream deps/addons that may have advanced since last release. edit recipes and psi4 `external/` code.
  * rebuild the psi4, psi4-rt, psi4-dev, psi4-docs stack on Linux and Mac

to do a rc release:
* bump the version and get it on master http://psicode.org/psi4manual/master/manage_git.html#how-to-bump-a-version
* build the explicit packages. don't rely upon usual build-latest but instead set branches directly at https://github.com/psi4/psi4meta/blob/master/conda-recipes/psi4-multiout/meta.yaml#L10 in psi4-multiout on L & M and https://github.com/psi4/psi4meta/blob/master/conda-recipes/psi4-docs/meta.yaml#L10 on L. start them building.

fill in codemeta.json
bump version https://github.com/loriab/psi4/blob/v13addlcontrib/codemeta.json#L9

*  add to the branch list https://github.com/psi4/psi4/blob/master/.travis.yml#L179
* bump copyright
  * tests/runtest.py
  * README.md
  * https://github.com/psi4/psi4/blob/fbb2ff444490bf6b43cb6e027637d8fd857adcee/doc/sphinxman/source/conf.py.in#L118
  * tests/psitest.pl

* (not done) main README updates

* tl;dr for version bump

```
# be on clean master up-to-date with upstream in both commits and tags
# * mind which version strings get "v" and which don't
# * if not fork, replace "upstream" with "origin"

>>> vi psi4/metadata.py
>>> git diff
diff --git a/psi4/metadata.py b/psi4/metadata.py
...
-__version__ = '1.3rc1'
-__version_long = '1.3rc1+5a7522a'
-__version_upcoming_annotated_v_tag = '1.3rc2'
+__version__ = '1.3rc2'
+__version_long = '1.3rc2+zzzzzzz'
+__version_upcoming_annotated_v_tag = '1.3rc3'

>>> git add psi4/metadata.py
>>> git commit -m "v1.3rc2"
[master bc8d7f5] v1.3rc2

>>> git log --oneline | head -1
bc8d7f5 v1.3rc2
>>> git tag -a v1.3rc2 bc8d7f5 -m "v1.3rc2"

>>> vi psi4/metadata.py 
>>> git diff
diff --git a/psi4/metadata.py b/psi4/metadata.py
...
-__version_long = '1.3rc2+zzzzzzz'
+__version_long = '1.3rc2+bc8d7f5'

>>> git add psi4/metadata.py
>>> git commit -m "Records tag for v1.3rc2"
[master 16dbd3e] Records tag for v1.3rc2

# goto GH:psi4/psi4 > Settings > Branches > master > Edit
#      https://github.com/psi4/psi4/settings/branch_protection_rules/424295
# uncheck "Include administrators" and Save changes

>>> git push upstream master
>>> git push upstream v1.3rc2

# re-engage "Include administrators" protections

```

build constructors by following instructions https://github.com/psi4/psi4meta/blob/master/conda-recipes/constructor-cutter-unified/README.md

* edit tag and hash
* clear package cache
* `python run.py`
* check for noarch (BAD)
* upload to vergil

generate the download site
* create new file akin to https://github.com/psi4/psicode-hugo-website/blob/master/content/installs/v13rc2.md mind the edition string `v13rc2` for this and future filenames
* copy and edit new https://github.com/psi4/psicode-hugo-website/blob/master/data/installs/v13rc2.yaml for menu and notes content
* enter `scripts/` dir and edit primarily https://github.com/psi4/psicode-hugo-website/blob/master/scripts/install-generator.py#L9 but also any other arrays or messages that should change.
* run the  `install-generator.py` in place. it will dump new files into `data/installs/` _subdirs_. be sure to `git add` them.
* installer page is ready. if wanted, edit the desination of the `Download` nav button https://github.com/psi4/psicode-hugo-website/blob/master/config.toml#L35
* pays to test
  * installer on Mac and Linux
  * that download button and `curl` downloading register on the download counters on vergil
  
  when publishing a new installer page, move the "latest" alias https://github.com/psi4/psicode-hugo-website/blob/master/content/installs/v13rc2.md to the new page so that "Downloads" on the naviagtion bar points there.
  
  in final pass before tag, check that external/ repos and commits match have been updated from conda recipes sources. also check versions with conda_build_config.yaml
check introduction.rst for any compiler and python min requirements to edit.

For final release:

### Repo/GH state

* follow tagging procedure
* before reengaing the "include admin" button, push a branch at the tag commit (not the records commit)

```
>>> git log --online | head -2
45315cb Records tag for v1.3
20e5c7e v1.3

>>> git checkout 20e5c7e
>>> git checkout -b 1.3.x
Switched to a new branch '1.3.x'
>>> git push upstream 1.3.x
```

* set up new branch as protected branch

### conda stage

* edit and switch to specific `git_tag` instead of `master` for psi4-multiout (L&M) and psi4-docs (L)
* in `ltrtver` in `conda_build_config.yaml`, make a new non-dev line (probably a ditto)

### conda stage (cont.)

* build the stack in the usual `***` configuration on L & M. should be (`psi4`, `psi4-rt`, `psi4-dev`) * python_versions for L & M. Also single `psi4-docs` from L.
* Poodle emits with `--label dev` so upload there to the subchannel directly. may need to clear out space.
* Go through each active conda package off https://anaconda.org/psi4/repo, find the most recent build set (L/M, active py versions) that psi/psi-rt/psi-dev is using and _add_ (not replace) the `main` label. this makes a `conda install psi4 -c psi4` get everything psi4 needs. for the moment `conda install psi4 -c psi4/label/dev` will get the same set, until package psi4-1.4a1.dev1 gets uploaded. may help to check versions and build versions against ltrtver in conda_build_config.yaml. this step takes a while.

### constructor stage

* move into constructor_cutter_unified. there's a rEADME there. 
  * edit pythons if necessary
  * edit release/hash/ltrtver of cookiecutter.json .
  * for non-rc's channel_tag should be empty string
  * leave git set to a `rc` since that has more details
  * copy cookiecutter.json to cookiecutter.json-vXXX
  * cookiecutter/{{.../construct.yaml rarely needs editing
  * do clear out .constructor so that everyting downloaded afresh
  * clean out build/
  * python run.py check for an py_ bad noarch pkgs
  * if fetching times out, may have to run run.py several times. clear out build/ in between. it's the fetching that takes a long time, not constucting
  * in the end, should have several constructors
```
>>> lh build/psi4conda-1.3-py3.*/*64.sh
-rwxr-xr-x. 1 psilocaluser psilocaluser 516M Feb 28 20:30 build/psi4conda-1.3-py3.6-linux-64/psi4conda-1.3-py36-Linux-x86_64.sh
-rwxr-xr-x. 1 psilocaluser psilocaluser 299M Feb 28 20:31 build/psi4conda-1.3-py3.6-osx-64/psi4conda-1.3-py36-MacOSX-x86_64.sh
-rwxr-xr-x. 1 psilocaluser psilocaluser 518M Feb 28 20:30 build/psi4conda-1.3-py3.7-linux-64/psi4conda-1.3-py37-Linux-x86_64.sh
-rwxr-xr-x. 1 psilocaluser psilocaluser 299M Feb 28 20:31 build/psi4conda-1.3-py3.7-osx-64/psi4conda-1.3-py37-MacOSX-x86_64.sh
```

* use command in readme to upload to vergil
* log in to vergil to make windows wsl symlinks

* worth downloading at least one L and M and installing it and running tests

generating download site directions above

shift alias to new installs/content file

NYI commit new files, PR, and upload site

On GH site "Draft a New Release"

grab the 1.3 manual by cp -pR master 1.3 on godaddy

NYI adjust the front-matter tags

reset
NYI new ltrt line
NYI psi4-multi and psi4-docs back to master
NYI turn off anom and go back to norm

edit RN and "publish" release. this establishes release data for GH api

close off RN issue

before stack build, consider max pinnings on deps, particularly any fast-moving deps (e.g., qcel) and whether they need version space to grow compatibly and grow incompatibly.

check in all release, constructor recipe changes on L/M. synchronize both to psi4meta

reset for normal operation

new ltrtver with new release.dev label
names back to master for psi4-multiout, psi4-docs
build string back to 0 if psi4-multiout needed multiple passes
poodle back to *** stack
crontab back to 2am "norm". comment out "anom"
new PR with edits to main README badges, py, etc

tweet
