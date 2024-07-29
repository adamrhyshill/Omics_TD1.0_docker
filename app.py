import dash
import dash_bootstrap_components as dbc

# meta_tags are required for the app layout to be mobile responsive
#app = dash.Dash(__name__, suppress_callback_exceptions=True,)
app = dash.Dash(__name__, suppress_callback_exceptions=True)
server = app.server

# sandstone
#pulse
# minty


#to run
# heroku login
# git init
# heroku git:remote -a halomonas-td10-omics
#halomonas-td10-omics
# git add .
# git commit -am "initial launch"
# git push heroku master
#

#to update the app, first close the app:
# heroku git:clone -a halomonas-td10-omics
# cd firstapp-omg
#next deploy your changes:
# git add .
# git commit -am "your update"
# source bin/activate


# heroku logs --tail to see error
# control c can get you out. 