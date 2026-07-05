# Releasing a new version of ShapeIt

This keeps the version shown in the About dialog, the GitHub tag/release, and
(once set up) the Zenodo archive all pointing at the same number, instead of
drifting independently the way "Version 1.1 (c) 2021" did.

Versioning follows [semver](https://semver.org): `MAJOR.MINOR.PATCH`
- **MAJOR**: incompatible settings-file format, or a change that alters
  analysis results for existing settings files
- **MINOR**: new features, backwards compatible
- **PATCH**: bug fixes only, no new features

## Steps

1. **Update `ShapeIt1.0/Include/ShapeVersion.h`**
   - Bump `SHAPEIT_VERSION` (e.g. `"1.2.0"` -> `"1.3.0"`)
   - Update `SHAPEIT_RELEASE_DATE` to today's date
   - Update `SHAPEIT_COPYRIGHT_YEARS` if the end year changed (e.g. `"2019-2026"` -> `"2019-2027"`)

2. **Commit the version bump on its own**, so it's easy to find later:
   ```bash
   git add ShapeIt1.0/Include/ShapeVersion.h
   git commit -m "Bump version to 1.3.0"
   git push
   ```

3. **Tag the commit, matching the version exactly** (with a `v` prefix, GitHub's convention):
   ```bash
   git tag -a v1.3.0 -m "v1.3.0"
   git push origin v1.3.0
   ```

4. **Create a GitHub Release from that tag** (Releases -> Draft a new release ->
   choose the tag you just pushed). Briefly describe what changed since the last
   release -- this becomes the changelog people see.
   - Anyone who has "Watch -> Custom -> Releases" set on the repo gets notified
     automatically at this step.

5. **If archiving to Zenodo for a DOI** (recommended once the software paper is
   submitted): if the repo is connected to Zenodo, creating the GitHub Release
   in step 4 triggers an automatic Zenodo archive of that exact tagged version,
   which is what you'd cite in the paper.

## Quick sanity check before tagging

Since the version is compiled into the code, rebuild and check
File -> About shows the version and date you just set, before pushing the tag.
