# CHANGELOG


## v0.4.0 (2025-04-10)

### Bug Fixes

- Update PyPI badge links to reflect the correct project name
  ([`187c3e9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/187c3e9f9f05c7af73ee5845128075ad4bbf26a5))

- Update README and .gitignore for consistency and clarity
  ([`1c6b0ed`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/1c6b0ed938beb3b021f7ea1193052f39431c6e21))

### Features

- Add pypi job for building and publishing packages
  ([`c50cc9f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c50cc9f6ecabcf44268e13f56339884b743524c5))

- Update project name and enhance metadata in pyproject.toml
  ([`3be44bb`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/3be44bbbe21ef9ab89b4a6cb45aa2cdae9ce3f83))


## v0.3.9 (2025-04-10)

### Bug Fixes

- Update pixi-pack-install-script version to v2
  ([`6489690`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/64896908c9eaf500b9091e14fd5160dedeadf725))


## v0.3.8 (2025-04-10)

### Bug Fixes

- Update checkout step to use release tag and bump pixi-pack-install-script version
  ([`4d481a7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/4d481a70fa77209663ff04bc968884b170b9976a))


## v0.3.7 (2025-04-10)

### Bug Fixes

- Add support for linux-aarch64 platform in release workflow
  ([`cb021ce`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/cb021cea133dff5eacd559562cdda53a3cc2f0f7))

- Update pixi-pack action version and adjust release workflow outputs
  ([`72c7fac`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/72c7fac91a9c874f8d89964b7f1c69f796a1406a))

### Documentation

- Add installation script command to README
  ([`71d7cc7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/71d7cc700dacfcdce47282181928dd58e200985c))


## v0.3.6 (2025-04-10)

### Bug Fixes

- Update release workflow to capture and use the release tag output
  ([`0738d37`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0738d37bd494d71cd9b6afe2f8351c9f2d5f7a35))


## v0.3.5 (2025-04-10)

### Bug Fixes

- Remove 'force' input from release workflow and adjust job conditions
  ([`a4ebebc`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a4ebebc4c5239d36c7dbb66f472b646c8980c688))

- Update pixi-pack action to v5 and modify build channels
  ([`0bb02a1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0bb02a1f93ece23de566799df057582d95fa83e3))


## v0.3.4 (2025-04-10)

### Bug Fixes

- Add 'force' input to release workflow to skip tests and benchmarks
  ([`c44bbd7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c44bbd7de66e53dbb8006e7ddf7142b2267afc92))

- Remove unused release commands and update semantic release configuration
  ([`8cdb251`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8cdb2517fd5ecbdfcb6e4e6528976e0275acec41))

- Update release job conditions to include 'force' input for immediate release
  ([`212ad06`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/212ad06e1d5ab72e1f2918b509383ba7f23c8356))

- Update SSH remote configuration and add repository URL for semantic release
  ([`6adbfb8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6adbfb89808fa87b43c1c77335df6e5d0dd85dba))


## v0.3.3 (2025-04-09)

### Bug Fixes

- Correct file extensions from .yml to .yaml in workflow configuration
  ([`af66615`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/af6661544c1a7ca6aded178eab87fcc2ddce1ca9))

- Correct file extensions in workflow configuration
  ([`00c9495`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/00c9495624b7542e2df5f0d5509ee3242e656b98))

- Correct permissions configuration in benchmark workflow
  ([`30430bf`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/30430bfeffa554c91e63c5823a1c8447196d510d))

- Enhance SSH configuration for deployment in release workflow
  ([`7a180e9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7a180e9ac564c4b1d4486f73d407c1a6d66cd037))

- Manual release workflow
  ([`2a6d406`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2a6d4065cd61ed90f973e7275411f4dfcefcd57b))

- Remove unnecessary checkout options in workflow configurations
  ([`a0e065d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/a0e065dfb26cc404e4fbef24fc55909d7e678c3c))

- Remove unnecessary git add and commit commands from release workflow
  ([`b2dbead`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/b2dbead5319eeed504a2f1e770fabd13eb50c93b))

- Remove unnecessary Git ref input from workflow configuration
  ([`50dba48`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/50dba48a7389ccdcdb8f46cfdadfb81a38f3943d))

- Remove unnecessary secrets from workflow configuration
  ([`054c62c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/054c62ccad2a8ab9668a99592bbe4de2b6eb05ce))

- Remove unnecessary secrets inheritance from benchmark workflow
  ([`f9d374a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f9d374a051aae6e96423e8512d4da84139ae2d33))

- Update job dependencies in release workflow
  ([`8cba440`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/8cba4408ccaa0a63d61cca20cf7215b2fead26ef))

- Update permissions and secrets configuration in workflows
  ([`2cf045c`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2cf045ce70a776152368ea3e818abddbdaf5160c))

- Update permissions to write for contents in workflows
  ([`f423c10`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/f423c1072ea33c1b390095c2f5849e6458d586bb))

- Update release command to prevent pushing during versioning
  ([`904903e`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/904903e7b78f27e4b818ee21cc77960b7431b8d7))

- Update release workflow to use semantic-release and rename jobs
  ([`6bb0728`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6bb0728be26b4e4c71307ec74e464a59440307ec))

### Chores

- Update commit message format for release step
  ([`ea4ced7`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/ea4ced70efb1c46269f952d0f6aff9bcae70751c))

### Continuous Integration

- Add build and release jobs to workflow
  ([`7dd8d9d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7dd8d9d24588f6d76eca823173c83a538f6975db))

### Refactoring

- Streamline CI workflows by modularizing tests and benchmarks
  ([`87ba582`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/87ba582965c7454293dd5f9a6d7c5377b5769465))


## v0.3.2 (2025-04-09)

### Bug Fixes

- Push changes
  ([`11816a9`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/11816a9e9971251c747a85053b12aab89a1cdf04))


## v0.3.1 (2025-04-09)

### Bug Fixes

- Commit release
  ([`4d494b8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/4d494b8998b27c675db78af3e87b8b72f548d39c))


## v0.3.0 (2025-04-09)

### Chores

- Remove benchmark workflow
  ([`c835035`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/c8350356cc3a70d7b1636173c29018d58a72dffa))

### Features

- Add benchmark to release workflow
  ([`0b9c3ff`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/0b9c3ff52b09dce9f4ff348245a250bc4e533579))

- Update release
  ([`6c3370d`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/6c3370d4116b345a7c5e53212b549a586819c0d6))


## v0.2.1 (2025-04-09)

### Bug Fixes

- Update release.yaml
  ([`36ffa42`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/36ffa42f5c47e92b98d5a48a2156f458f3c5e4c1))


## v0.2.0 (2025-04-09)

### Features

- Add ssh key
  ([`814c2d8`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/814c2d8c80aee7006c07b2e243cfcc8d704eed18))

- Add token to release workflow
  ([`2e3b879`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/2e3b8796eca56869580b45d8716c774e93f88ac9))


## v0.1.0 (2025-04-09)

### Bug Fixes

- Replace semantic-release action
  ([`005dab1`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/005dab184ff80719f55ca51b40224a90a2ad9c06))

### Features

- Add manual release workflow
  ([`d01c77a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/d01c77a97b1f772a605b8237f2efb5d34774113d))


## v0.1.0-rc.1 (2025-04-09)

### Bug Fixes

- Rename mkdocs.yml to mkdocs.yaml
  ([`fbbf84f`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/fbbf84fb0a811bec5dc61f2fa8b94c7d519b8191))

- Rename workflow file for PRs
  ([`7ff1490`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/7ff1490310248fb99f3643e6873c63f04b7c9cd4))

### Features

- Add release workflow
  ([`dba2b3a`](https://github.com/centre-pathogen-genomics/snippy-ng/commit/dba2b3a790a67a68e4019cc7c6b61610575a8bda))
