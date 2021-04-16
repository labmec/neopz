/**
 * @file
 * @brief Contains declaration of TPZGuiInterface class.
 */
//$Id: TPZGuiInterface.h,v 1.2 2011-04-05 19:32:54 calle Exp $

#ifndef TPZGuiInterfaceH
#define TPZGuiInterfaceH

#include "pzerror.h"
#include "TPZSavable.h"
/**
 * @ingroup common
 * @brief This class implements a very simple interface from PZ kernel to GUI. Module: \ref common "Common".
 * @author Tiago Forti
 * @since 2010, March 30
 */
/** 
 * It implements for instance the capability of cancel the execution
 * The GUI must define a derived class which reimplements the messages and update methods
 * for better messages.
 */
class TPZGuiInterface : public TPZSavable {
	
protected:
	
	/** Flag indicating if the thread is alive or if its canceling was requested */
	bool fCanceled;
	
	/** Message attribute for UpdateGUI method */
	std::string fMessage;
	
	/** Progress bar attributes for UpdateGUI method */
	int fProgressBarPos, fProgressBarMaxPos, fProgressBarMinPos;
	
public:
	
	/** Default constructor */
	TPZGuiInterface();
	
	/** Default destructor */
	virtual ~TPZGuiInterface();
    
    TPZGuiInterface(const TPZGuiInterface &gui) : fCanceled(gui.fCanceled), fMessage(gui.fMessage), fProgressBarPos(gui.fProgressBarPos),
        fProgressBarMaxPos(gui.fProgressBarMaxPos),fProgressBarMinPos(gui.fProgressBarMinPos)
    {
        DebugStop();
    }
    
    TPZGuiInterface &operator=(const TPZGuiInterface &gui)
    {
        fCanceled = gui.fCanceled;
        fMessage = gui.fMessage;
        fProgressBarPos = gui.fProgressBarPos;
        fProgressBarMaxPos = gui.fProgressBarMaxPos;
        fProgressBarMinPos = gui.fProgressBarMinPos;
        DebugStop();
        return *this;
    }
	    int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;
    
	/** Updates the GUI with start messages \n
	 * This method must be reimplemented in derived classes for better messages
	 */
	virtual void Start();
	
	/** Updates the GUI \n
	 * This method must be reimplemented in derived classes for better messages
	 */
	virtual void UpdateCaption();
	
	/** Updates the GUI with ending messages \n
	 * This method must be reimplemented in derived classes for better messages
	 */
	virtual void End();
	
	/** Show an error message. \n
	 * This method must be reimplemented in derived classes for better messages
	 */
	virtual void ShowErrorMessage(std::string message);
	
	/** Defines the process has been canceled */
	void SetKilled();
	
	/** Returns whether the process was canceled */
	bool AmIKilled();
	
	/** Message attribute for UpdateGUI method */
	std::string Message();    
    /** Change the message of the object */
    void SetMessage(const std::string &message);
	
	/** Progress bar attributes for UpdateGUI method */
	int &ProgressBarPos(){
		return fProgressBarPos;
	}
	
	/** Progress bar attributes for UpdateGUI method */
	int &ProgressBarMaxPos(){
		return fProgressBarMaxPos;
	}
	
    /** Progress bar attributes for UpdateGUI method */
	int &ProgressBarMinPos(){
		return fProgressBarMinPos;
	}
	
};

//---------------------------------------------------------------------------
#endif
