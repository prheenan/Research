Attribute VB_Name = "NewMacros"
Sub wrap_table()
Attribute wrap_table.VB_ProcData.VB_Invoke_Func = "Normal.NewMacros.Macro3"
'
' wrap_table Macro
'
'
    ' Make an index array; only want to align the first N tables
    '
    Selection.WholeStory
    Dim i As Integer
    i = 0
    ' Loop through each table
    '
    Dim t As Table
    For Each t In ActiveDocument.Tables
        i = i + 1
        If (i = 1 Or i = 2 Or i = 3) Then
            t.Rows.WrapAroundText = True
            t.Rows.HorizontalPosition = wdTableLeft
            t.Rows.VerticalPosition = InchesToPoints(1)
            t.Rows.RelativeVerticalPosition = wdRelativeVerticalPositionPage
            t.Rows.AllowOverlap = False
        End If
    Next
    'Biophysical journal wants...
    ' (1) justified paragraphs
    Selection.ParagraphFormat.Alignment = wdAlignParagraphJustify
    ' (2) Time New Roman Font
    Selection.Font.Name = "Times New Roman"
End Sub
