def version_scheme(version):
    return "1.0.0b0"


def local_scheme(version):
    return f"git.{version.node[:7]}"
