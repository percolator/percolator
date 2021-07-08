# Code Review Process

To begin with, it is recommended that you develop your code in your
own fork of the Github repository, rather than commit directly to the
main repository. Github makes it easy to submit code reviews using
this workflow. This will also allow you to validate your change for
the various builds before submitting your changes for review.

(The alternative is to commit your changes in the main repository but
on a separate branch. You can then request a code review before
merging your branch to the master branch. This may be more appropriate
in some cases.)

## Submitting a Code Review

Once your code has been pushed to your own repository's master branch,
and the build action verifies that the code changes build correctly,
you can submit it for review via the following steps:

1. Visit your repository's web page on `github.com`. At the top of the
directory, there will be a line indicating that your repository is now
"ahead of" the original repository by one or more commits. There will
be a button labeled "Contribute". Use this to pull down a button
labeled "Open pull request", and select it.

2. On the resulting page, you will be shown a diff of your changes.
Scan the page to verify that it looks correct, and doesn't include
anything unexpected (such as overlooked debugging lines or extraneous
whitespace). Assuming all is well, select the "View pull request"
button to proceed.

3. On the far right of the resulting page, there will be a section
labeled "Reviewers". Select the gear icon to pull down a list of
people that can review your changes. Select an appropriate reviewer.

4. While on this page, check that the title of the pull request is a
one-line summary the change to the reviewer. By default this title
will be taken from your last commit, which may not be a useful summary
of the overall change. If a fuller description is called for, use the
description box to provide details.

5. When everything is ready, use the "Create pull request" button to
submit your change for review.

## Reviewing a pull request

When a code reivew has been requested from you, you will receive an
email notification via Github. The code change will also be in the
repository's list of open pull requests. Follow the link from either
of these places to view the change. You can explicitly select the
"Files changed" link from the summary page to bring up the complete
diff for this pull request.

Read through the diff. If and when you see a change that warrants a
comment from you, hover over the relevant line in the green, changed
sectionk, and a small button with a plus sign will appear at the start
of the line. Select this button to bring up an edit control, in which
you can make your comment. Use the button labeled "Start a review" to
add the comment to your review. (The next comment that you compose
will get a button labeled "Add review comment".)

When you have finished examining the change, there will be button back
up top labeled "Finish your review". This will drop down a control
that will allow you to make an overall comment, if necessary. You can
then choose to "Comment", "Approve", or "Request changes". Usually you
will choose "Approve", unless your comments include an observation
that is uncertain or non-trivial to address. ("Comment" is the proper
choice if you are not the main reviewer for this pull request.) You
can then select "Submit review" to have your comments sent to the
requester.

If you found no issues in your review, and the change is a
straightforward one, you may choose to merge the pull request
yourself. However, it is more common to let the original requester do
the actual merge, as they are the one who is ultimately responsible
for the change.
