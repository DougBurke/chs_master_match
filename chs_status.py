"""
Indicate ensemble status.

Perhaps this should be handled by an object rather than simple
labels?

Note that QA cases are actually "qa-<reason>", and ideally
we'd identify that here. For now there is the is_qa function.

"""

__all__ = ("COMPLETED", "DELETE", "DONE", "FINALIZE", "QA",
           "REVIEW", "TODO", "UNKNOWN", "is_qa")


COMPLETED = "completed"
DELETE = "delete"
DONE = "done"
FINALIZE = "finalize"
QA = "qa"
REVIEW = "review"
TODO = "todo"
UNKNOWN = "unknown"


def is_qa(status):
    """Is the status a QA status?

    Parameters
    ----------
    status : str
        The status string.

    Returns
    -------
    flag : bool
        True if status begins with 'qa-'.

    """

    return status.startswith('qa-')
