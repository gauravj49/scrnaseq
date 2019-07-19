# Initialize a repository from existing code (Ex. TRN)
# Add the files in your new local repository. This stages them for the first commit.
git init
# Initialized empty Git repository in /home/rad/users/gaurav/projects/trn/.git/

# Set config values
git config user.name  "Gaurav Jain"
git config user.email "gauravj49@gmail.com"
touch README.md

# Check the config list
git config --list

# Get the status of the project and repository
git status

# Ignore files that should not go into the repository
emacs -nw .gitignore
# annotation
# docs
# input
# output
# .gitignore
# 00_git_setup_usage.sh

# Adds the files in the local repository and stages them for commit. To unstage a file, use 'git reset HEAD YOUR-FILE'. Commit the files that you've staged in your local repository.
git add -A

# Remove file from staging area
git reset <optional filenames>

# Make the first commit
git commit -m "Initial commit of TRN (transcriptional regulatory network) identification code"

# Add to github
# Create a project on github and do not initialize it. We have already done that here.
# Source: https://stackoverflow.com/questions/12799719/how-to-upload-a-project-to-github
# So far, the above steps is what you would do even if you were not using github. They are the normal steps to start a git repository. Remember that git is distributed (decentralized), means you don't need to have a "central server" (or even a network connection), to use git. Now you want to push the changes to your git repository hosted with github. To you this by telling git to add a remote location, and you do that with this command: git remote add origin https://github.com/yourusername/your-repo-name.git
# curl -u 'USER' https://api.github.com/user/repos -d '{"name":"REPO"}'

curl -u 'gauravj49' https://api.github.com/user/repos -d '{"name":"TRN"}'
git remote add origin https://github.com/gauravj49/TRN.git

# Once you have done that, git now knows about your remote repository. You can then tell it to push (which is "upload") your commited files:
git push -u origin master