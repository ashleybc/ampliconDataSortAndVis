import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import colorcet as cc
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.markers as mmarkers
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import scipy as sp
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as plticker


#my functions
def topCruiseRelPicks(df,cruises,noPicks):
#This function returns a list of unique taxa from the lists of relative top taxa per cruise
#for example: The most relatively abundant taxa in all GL3 samples is some list of 20 taxa, and then
#there is a similar list for GL4 samples, and then there are 20 unique taxa among them, which is then what we plot.
#must go with weighted/rel abundances because top picks by total reads will be biased towards samples with more reads
#rather than capturing what is the most relatively important at all depths.

#i found that setting noPicks to 15 for everything works well (~40 genera, 30 orders, 25 classes)
    allTaxa=[]
    for c in cruises:
        #print("cruise",c)
        
        colList=[]
    
        for col in df:
            if c in col:
                #print(col)
                colList.append(col)
                #print(colList)
        #print("\n")
        sub=df[colList]
        subtaxasum=sub.sum(axis=1)
        #print(sub)
        #print("\n")
        #print(subtaxasum)
        #print("\n")
        sortedsums=subtaxasum.sort_values(ascending=False)[:noPicks]
        #print(sortedsums)
        sortedsumTaxa=list(sortedsums.index)
        
        for t in sortedsumTaxa:
            allTaxa.append(t)
            
    finList= list(np.unique(allTaxa))
    
    return finList

def plt_names(pltdf,lev):

    plotNames=[]
    for col in pltdf.columns:
    
        if "D_0__Archaea" in col:
            arcTag=" (A)"
            
        else:
            arcTag=""
            
        pltnamestart=col.find(";D_"+str(lev-1)+"__")
        levUpstart=col.find(";D_"+str(lev-2)+"__")
        
        if not pltnamestart==-1 and col != "other":
            pltname=col[pltnamestart+len(";D_"+str(lev-1)+"__"):]
            levUpName=col[levUpstart+len(";D_"+str(lev-2)+"__"):pltnamestart]
            #print("start pltname", pltname)
            #print("prior lev",levUpName)
            
            if levUpName.lower().startswith("clade"):#case sensitive, so make sure everything is one case
                levUpName=col[col.find(";D_"+str(lev-3)+"__")+len(";D_"+str(lev-3)+"__"):levUpstart]+" "+levUpName
                #print("new prior lev",levUpName)
            for i in range(lev)[::-1]:
            
                if "uncultured" in pltname or "bacterium" in pltname or "ambiguous" in pltname.lower():
                    pltnamestart=col.find(";D_"+str(i-1)+"__")
                    pltnameend=col.find(";D_"+str(i)+"__")
                    pltname="Unassigned "+col[pltnamestart+len(";D_"+str(i-1)+"__"):pltnameend]
                else:
                    break
                
        elif col !="Other":
            firstDiv=col.rfind("D_")
            pltname="Unassigned "+col[firstDiv+len("D_")+3:]
            levUpName=""
        else:
            pltname=col
            levUpName=""
            
        if lev !=2:
            if "Unassigned" in pltname and not "clade"in levUpName.lower() or pltname=="Other":
                pltname=pltname.rstrip(";__")
            else:
                pltname=(levUpName+";"+pltname).rstrip(";__")
        else:
            pltname=pltname.rstrip(";__")
        
        plotNames.append((pltname+arcTag).replace("_"," "))
        
        uniquePlotNames=list(np.unique(plotNames))
        
        #check for/fix non-unique names
        for u in uniquePlotNames:
            if plotNames.count(u)>1:
                #print(u)
                ndices = [i for i, x in enumerate(plotNames) if x == u]
                #print(ndices)

                for count,n in enumerate(ndices):
                    plotNames[n]=u+" "+str(count+1)

    return plotNames

def plot_pre_proc(df,taxa,lev,md,norm=False,fullchangelab=True):

    plotdf=df.loc[taxa].T #transpose dataframe so taxa are columns and trap ID is index
    plotdfsums=plotdf.sum(axis=1)
    
    if norm==True:
        plotdf.name="top"+df.name+"norm"
    else:
        plotdf.name="top"+df.name

    if norm==True:
        for ind in plotdf.index:
            #print(plotdf.loc[ind])
            #print(plotdfsums.loc[ind])
            plotdf.loc[ind]=plotdf.loc[ind].apply(lambda x:x/plotdfsums.loc[ind])

    
    else:
        plotdf["Other"]=plotdfsums.apply(lambda x: 1.00000-x) #make category "other for whatever doesn't make it into top taxa
        #print(plotdf["Other"])
    
    plotdf.columns=plt_names(plotdf,lev)
    
    return plotdf


def stacked_bar_view(plotdf,lev,norm=False):

    newIDs=[]
    for ind in plotdf.index:
        newIDs.append(metadata.loc[ind]["Month"].replace("_"," ")+" "+str(int(metadata.loc[ind]["Year"]))+", "+str(int(metadata.loc[ind]["Depth"]))+" m")
    
    plotdf.index=newIDs
    plotdf.sort_index(inplace=True)

    if norm==True:
        normStr="norm"
    else:
        normStr=""

    fig = plt.figure()
    ax1 = plt.subplot(111)
    
    plotdf.plot.bar(stacked=True,ax=ax1,cmap=cc.cm.glasbey_light,edgecolor="black")
    plt.subplots_adjust(right=0.6)

    plt.ylim(0,1.0)
    
    if norm==True:
        ax1.set_ylabel("Adjusted Relative Abundance")
    else:
        ax1.set_ylabel("Relative Abundance")
    # Put legend to the right of the current axis
    ax1.legend(fontsize=9,loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.savefig(root+"/pythoncode/trapseqorg1out/top"+str(len(plotdf.columns))+"L"+str(lev)+normStr+".pdf",bbox_inches="tight")
    plt.close()

def bubble_pre_proc(plotdf,lev,cruises):

    plotdfcopy=plotdf.copy()
    
    cruiseOrder=[] #create column to do .groupby to retreive data by cruise
    for ind in plotdfcopy.index:
        #print(ind)
        for c in cruises:
            if ind.startswith(c):
                cruiseOrder.append(c)
                break
    
    #print("cruise order",cruiseOrder)
    plotdfcopy["cruise"]=cruiseOrder
    plotdfcopy["Depth"]=metadata.loc[plotdfcopy.index]["Depth"]
    plotdfcopy=plotdfcopy.astype({"Depth":int})
    
    plotdfcopy.name=plotdf.name
    
    uniqueDepths=list(np.unique(plotdfcopy["Depth"]))

    subs=[]

    for count,g in enumerate(plotdfcopy.groupby("cruise")):
        #print(count)
        
        temp=g[1] #data frame is second position in tuple of group,df
        #print(temp)
        dps=list(temp["Depth"]) #depths in dataframe group
        #print(dps)
        missing=[i for i in uniqueDepths if i not in dps] #determine which depths of all depths are not there
        #print(missing)

        fakeRows=[]
        fakeInds=[]
        for m in missing: #for each missing depth
            fakeInd=g[0]+"-"+str(int(m)) #create fake index and row
            fakeDepth=m
            fakeCruise=g[0]
            fakeFiller=[np.nan]*len(temp.columns)

            fakeFiller[-1]=fakeDepth
            fakeFiller[-2]=fakeCruise
            
            fakeRows.append(fakeFiller)
            fakeInds.append(fakeInd)

        addDf=pd.DataFrame(fakeRows, columns=list(temp.columns),index=fakeInds) #make dataframe of fake data

        filledDf=temp.append(addDf) #append to existing
        filledDf.sort_values("Depth",inplace=True) #sort by depth shallow-->deep
        filledDf=filledDf.astype({"Depth":str}) #change to string to make labels=exact depths
        filledDf.drop("cruise",axis=1,inplace=True)
        
        if plotdfcopy.name.endswith("norm"):
            filledDf.name=g[0]+"L"+str(lev)+"norm"
        else:
            filledDf.name=g[0]+"L"+str(lev)
            
        subs.append(filledDf)
        
    return subs

def customLegend(legColList,fillStList,legLabList,markerList):

    hList=[]

    for count, s in enumerate(legColList):

        h=mlines.Line2D([], [], color=s, marker=markerList[count], fillstyle=fillStList[count],linestyle="None",markersize=13, label=legLabList[count],markeredgecolor="None",alpha=0.5)
        hList.append(h)
    
    return hList

def catbincolors(fr,binfr,cat,colorlist):

    binseps = eval(binfr.loc[cat].bin_separators)

    fcList=[]
        
    for d in fr[cat]:
        print(d)
        for ind in range(len(binseps)+1):
            if ind==0:
                #print("<"+str(binseps[ind]))
                if int(d) < binseps[ind]:
                    #print("True")
                    fcList.append(colorlist[ind])
                    break
                        
            elif ind==len(binseps):
                #print(">="+str(binseps[-1]))
                if int(d) >= binseps[-1]:
                    fcList.append(colorlist[ind])
                    break
            else:
                #print("<"+str(binseps[ind])+">="+str(binseps[ind-1]))
                if binseps[ind] > int(d) >= binseps[ind-1]:
                    fcList.append(colorlist[ind])
                    break

    return fcList


def bubble_plot(subsList,lev,subplLabs,cat,catbins=pd.DataFrame([]),catcolors=[],labSuffix="",monthreord=False):

    if subsList[0].name.endswith("norm"):
        normStr="norm"
    else:
        normStr=""
    
    
    if monthreord==True:
        reord=["GL1","GL4","GL2","GL3"]
        
        currentOrd=[]
        
        for s in subsList:
            currentOrd.append(s.name[:s.name.rfind("L")])
        
        reorderedSubs=[]
        for n in reord:
            ind=currentOrd.index(n) #get the index of term in desired order list in the old one
            toappend=subsList[ind] #get df by index
            reorderedSubs.append(toappend)
            
        #print("reordered: ")
        #for r in reorderedSubs:
            #print(r.name)
            
        subsList=reorderedSubs
    
    plt.rc("text", usetex=True)
    #begin the figure, set the number of subplots equal to the number of dataframes
    fig,axes = plt.subplots(nrows=1, ncols=len(subsList), figsize=(15, 9), sharey=True)

    #for each dataframe, make a copy so that the original dataframe is not modified
    for count,fr in enumerate(subsList):
        #print(fr.name)
        #print(fr[cat])
        
        if catbins.empty==False:
            fcList=catbincolors(fr,catbins,cat,catcolors)
        else:
            fcList=["blue"]*len(fr[cat])

        for colcount,col in enumerate(fr.columns):
            if col !=cat and col !="Other" and col != "other":
                axes[count].scatter(x=fr[cat],y=[col]*fr.index.size,s=fr[col]*3000, facecolor=fcList,marker=mmarkers.MarkerStyle("o", fillstyle="full"),alpha=0.5,edgecolor=fcList,linewidths=0.8)
        
        axes[count].grid(linestyle=":",axis="both")
        axes[count].set_axisbelow(True)
        axes[count].set_xlim(-1,len(fr[cat]))
        
        axes[count].xaxis.set_label_position("top")
        axes[count].set_xlabel(subplLabs[count],fontsize=12)
        
        for tick in axes[count].get_xticklabels():
            tick.set_rotation(70)
            
    #colors legend
    if catbins.empty==False:
        l1handList=customLegend(catcolors,["full"]*len(catcolors),list(myPosColorDict.keys()),["o"]*len(catcolors))

        plt.subplots_adjust(wspace=0.1)

        #must use fig.legend rather than plt.legend to place legend in figure space outside of subplots
        legend1=fig.legend(handles=l1handList,labels=list(myPosColorDict.keys()),fontsize=12,ncol=3,bbox_to_anchor=[0.2,0.001],loc="lower left")
        fig.add_artist(legend1)
    
    plt.text(0.045, 0.08, cat+labSuffix, fontsize=12, transform=plt.gcf().transFigure)
    plt.savefig(root+"/pythoncode/trapseqorg1out/top"+str(len(subsList[0].columns))+"L"+str(lev)+normStr+"bubble.pdf",bbox_inches="tight")
    plt.close()


#other user inputs
taxdict={0:'Domain',1:'Phylum',2:'Class',3:'Order',4:'Family',5:'Genus',6:'Species'}
myLevs=[2,3,4,5,6]
cruises=["GL1","GL2","GL3","GL4"]
dropCols=["GEN-DONOR","MOCK-EVEN", "MOCK-STAG","tr-blk"]
#file import
root="/Users/ashley/Documents/Research/Gordon'sLab/FGL/amplicon"

metadata=pd.read_csv(root+"/Qiime2_analyzed_FGL_traps/metadata.tsv",delimiter="\t",index_col="SampleID")
metadata.drop(dropCols,axis=0,inplace=True)

legendBins=pd.read_csv(root+"/Qiime2_analyzed_FGL_traps/legendBins.txt", sep="\t",index_col=0)

myPosColorDict={"mixo.+redoxcline":"green","plate":"purple","mid+deep monim.":"red"}

def runner():
    dfs=[]
    topLevBarPicks=[]
    for lev in myLevs:
        #print("input level",lev)
        df=pd.read_csv(root+"/Qiime2_analyzed_FGL_traps/rel_OTUtables/rel-phyla-table_L"+str(lev)+".tsv", delimiter="\t",skiprows=[0],index_col="#OTU ID")
        df.drop(dropCols,axis=1,inplace=True)
        df.name="allrelL"+str(lev)
        dfs.append(df)

        if lev==2:
            #print(lev)
            barPicks=topCruiseRelPicks(df,cruises,10)

        elif lev==6:
            barPicks=topCruiseRelPicks(df,cruises,12)

        else:
            barPicks=topCruiseRelPicks(df,cruises,15)

        topLevBarPicks.append((lev,barPicks))

    topLevBarTaxaDict=dict(topLevBarPicks)

    procdfs=[]
    for d in dfs:
        lev=d.name[-1]
        procdf=plot_pre_proc(d,topLevBarTaxaDict[int(lev)],int(lev),metadata,norm=False)
        procdfnorm=plot_pre_proc(d,topLevBarTaxaDict[int(lev)],int(lev),metadata,norm=True)
        
        procdfbubble=plot_pre_proc(d,topLevBarTaxaDict[int(lev)],int(lev),metadata,norm=False,fullchangelab=False)
        procdfbubbleNorm=plot_pre_proc(d,topLevBarTaxaDict[int(lev)],int(lev),metadata,norm=True,fullchangelab=False)
        
        stacked_bar_view(procdf,lev,norm=False)
        
        bubbledfsubs=bubble_pre_proc(procdfbubble,lev,cruises)
        bubbledfsubsNorm=bubble_pre_proc(procdfbubbleNorm,lev,cruises)
        
        bubble_plot(bubbledfsubs,lev,["early Jul. 2016","late Jul. 2017","early Aug. 2016","early Oct. 2016"],"Depth",legendBins,list(myPosColorDict.values())," (m)",monthreord=True)
        bubble_plot(bubbledfsubsNorm,lev,["early Jul. 2016","late Jul. 2017","early Aug. 2016","early Oct. 2016"],"Depth",legendBins,list(myPosColorDict.values())," (m)",monthreord=True)

