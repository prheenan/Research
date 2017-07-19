Attribute VB_Name = "NewMacros"
Sub spell_check()
Attribute spell_check.VB_ProcData.VB_Invoke_Func = "Normal.NewMacros.spell_check"
'
' spell_check Macro: spell checks the document
'
'
    Selection.WholeStory
    Selection.LanguageID = wdEnglishUS
    Selection.NoProofing = False
End Sub
