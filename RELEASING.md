# Releasing PyCoast

1. checkout master
2. pull from repo
3. run the unittests
4. run `loghub` and update the `CHANGELOG.md` file:

```
loghub pytroll/pycoast -u <username> -st v<last-tag> -plg bug "Bugs fixed" -plg enhancement "Features added" -plg documentation "Documentation changes"
```

Don't forget to commit!

5. Create a tag with the new version number, starting with a 'v', eg:

```
git tag v0.22.45
```

See [semver.org](http://semver.org/) on how to write a version number.



6. push changes to github `git push --follow-tags`

Make sure the new tag has been pushed.

7. Verify travis tests passed and deployed sdist to PyPI
