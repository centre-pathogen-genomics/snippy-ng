from glasscandle import Watcher
from glasscandle.notifications import slack_notifier


slack_notify = slack_notifier()

watch = Watcher(
    "dependencies/versions.json",
    on_change=slack_notify
)

watch.conda("bcbio-gff", channel="bioconda")
watch.conda("click", channel="conda-forge")
watch.conda("packaging", channel="conda-forge")
watch.conda("pydantic", channel="conda-forge")
watch.conda("biopython", channel="conda-forge")
watch.conda("samtools", channel="bioconda")
watch.conda("samclip", channel="bioconda")
watch.conda("bwa", channel="bioconda")
watch.conda("freebayes", channel="bioconda")
watch.conda("bcftools", channel="bioconda")
watch.conda("vt", channel="bioconda")


if __name__ == "__main__":
    updated = watch.run()
