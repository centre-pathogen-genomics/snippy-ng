import tomllib

from glasscandle import Watcher
from glasscandle.notifications import slack_notifier


slack_notify = slack_notifier()
watch = Watcher("dependencies/versions.json",on_change=slack_notify)
# watch = Watcher("dependencies/versions.json")

# dynamically add dependencies from pyproject.toml
with open('pyproject.toml', 'rb') as f:
    data = tomllib.load(f)

for dep, details in data['tool']['pixi']['dependencies'].items():
    version = details.get('version')
    version = None if version == "*" else version # * is not a valid version constraint
    watch.conda(dep, version=version, channels=[details['channel']])

if __name__ == "__main__":
    updated = watch.run()
