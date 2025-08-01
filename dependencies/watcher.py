import tomllib

from glasscandle import Watcher
from glasscandle.notifications import slack_notifier


slack_notify = slack_notifier()

watch = Watcher(
    "dependencies/versions.json",
    on_change=slack_notify
)

# dynamically add dependencies from pyproject.toml
with open('pyproject.toml', 'rb') as f:
    data = tomllib.load(f)

for dep in data['tool']['pixi']['dependencies']: 
    watch.conda(dep)

if __name__ == "__main__":
    updated = watch.run()
