
import os, sys, re, time

print('\n')

print('This script will initialize a project at {0}. It will create the following folders:\n'.format(os.getcwd()))

print("- A remote folder, which is a symbolic link to a folder where you store the data for this project, ussually located in your Dropbox or Google Drive. This folder will be ignored by git.\n")

print("- A local folder where you can locally store data you don't want to have in your remote folder (e.g. large files that are used in intermediate steps). You should programmatically create the necessary subdirectories and files if you want your colaborators to be able to reproduce your results with the shared code. This folder will be ignored by git.\n")

print("- A script folder where you can put the scripts you share with your colaborators. This folder will be tracked by git.\n")

print("- A tex folder where you can put the tex files you share with your colaborators. This folder will be tracked by git.\n")

print("You can costumize the .gitignore file. As it is now, it will ignore the directories mentioned above, the standard python files you don't want to track, and the usual latex files you don't want to track.\n")

print("You will be ask now for the path to your remote folder. If a folder with the same name already exists, the script will skipped it and warns the user.\n")

time.sleep(1)

if os.path.exists('remote'):
    print('A remote directory already exists. Skipping')
else:
    while True:
        target_dir = input('Please provide the absolute path to the remote directory in your system: ')
        target_dir = target_dir.replace("'",'').replace('"','')
        if not os.path.exists(target_dir):
            print('{0} was not found in your system. You may include the path by dragging the directory into the terminal, but you have to delate any extra spaces added at the end.'.format(target_dir))
            continue
        else:
            os.symlink(target_dir,os.path.join('remote'))
            break

if os.path.exists('local'):
    print('A local directory already exists. Skipping')
else:
    os.mkdir(os.path.join('local'))

if os.path.exists('scripts'):
    print('A scripts directory already exists. Skipping')
else:
    os.mkdir(os.path.join('scripts'))

if os.path.exists('tex'):
    print('A tex directory already exists. Skipping')
else:
    os.mkdir(os.path.join('tex'))
    os.mkdir(os.path.join('tex','tables'))
    os.mkdir(os.path.join('tex','figures'))
