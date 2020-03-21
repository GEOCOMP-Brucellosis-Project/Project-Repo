# Readme!

Maybe this can be a message board of sorts



## Setting up an environment
To run the code you will need all of the packages installed. One way to do this would be to try running the code and install packages one-by-one 
based on the errors that python gives you. This will probably work, however, sometimes it is better to install certain packages before others. 
To install packages with dependencies in mind, you can install them all at once.

I made a list of the packages that I am currently using and uploaded it. To create an env identical to the one that I have for running the code
place spec-file.txt in your working directory and run the following in command prompt:
conda create --name "whatever you want to call your env" --file spec-file.txt

Then follow instructions here to activate the environment/etc. 
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#building-identical-conda-environments


When we're writing code lets just keep each other updated if we use a new package in the project. 



## EarthEngine Credentials

When you call ee.Initialize() it searches for a credentials file on your computer. If you use mine, then I think you should be able to access my assets (like the iran shapefile.)
Otherwise we might need to use different code for this part, since earth engine only plays nicely with "assets" hosted on earth engine. You can try putting the credentials file
in the spot shown on the screenshot. 